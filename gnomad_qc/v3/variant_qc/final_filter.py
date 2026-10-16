# noqa: D100

import argparse
import logging
import sys
from typing import Optional

import os
import hail as hl

from gnomad_qc.pipeline_config import BASE, WRITE_DIR, init_hail

init_hail("final_filter")

from gnomad_qc.v3.resources.filtering import add_filters_expr
from gnomad.resources.grch38.reference_data import telomeres_and_centromeres
from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists
#from gnomad.utils.filtering import add_filters_expr
from gnomad.utils.slack import slack_notifications
from gnomad.variant_qc.pipeline import INBREEDING_COEFF_HARD_CUTOFF

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.annotations import get_freq, get_info, get_vqsr_filters
from gnomad_qc.v3.resources.release import release_sites
from gnomad_qc.v3.resources.variant_qc import final_filter, get_score_bins

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variant_qc_filtering")
logger.setLevel(logging.INFO)


def generate_final_filter_ht(
    ht: hl.Table,
    model_name: str,
    score_name: str,
    inbreeding_coeff_cutoff: float = INBREEDING_COEFF_HARD_CUTOFF,
    snp_bin_cutoff: int = None,
    indel_bin_cutoff: int = None,
    snp_score_cutoff: float = None,
    indel_score_cutoff: float = None,
    aggregated_bin_ht: Optional[hl.Table] = None,
    bin_id: Optional[str] = None,
    vqsr_ht: Optional[hl.Table] = None,
) -> hl.Table:
    """
    Prepare finalized filtering model given a filtering HT from `rf.apply_rf_model` or VQSR.
    """

    # --- Step 1: ensure InbreedingCoeff exists ---
    ht = ht.annotate(InbreedingCoeff_final=hl.coalesce(ht.InbreedingCoeff, ht.InbreedingCoeff_1))

    # --- Step 2: compute filter expressions directly on this table ---
    ht = ht.annotate(
        ac0_filter_expr=hl.or_else(ht.freq_ac0 == 0, False),
        monoallelic_flag_expr=hl.or_else((ht.freq_monoallelic_af == 1) | (ht.freq_monoallelic_af == 0), False),
        ts_ac_filter_expr=hl.or_else(ht.freq_transmitted_singleton_ac == 1, False),
        inbreeding_coeff_flag=ht.InbreedingCoeff_final < inbreeding_coeff_cutoff,
        transmitted_singleton_from_ts=ht.freq_transmitted_singleton_ac == 1
    )

    # --- Step 3: Validations for SNP/indel cutoffs --- 
    if snp_bin_cutoff is not None and snp_score_cutoff is not None: 
        raise ValueError("snp_bin_cutoff and snp_score_cutoff are mutually exclusive. Supply only one SNP cutoff.") 
    if indel_bin_cutoff is not None and indel_score_cutoff is not None: 
        raise ValueError("indel_bin_cutoff and indel_score_cutoff are mutually exclusive. Supply only one indel cutoff.") 
    if snp_bin_cutoff is None and snp_score_cutoff is None: 
        raise ValueError("Supply one of snp_bin_cutoff or snp_score_cutoff.") 
    if indel_bin_cutoff is None and indel_score_cutoff is None: 
        raise ValueError("Supply one of indel_bin_cutoff or indel_score_cutoff.") 
    if (snp_bin_cutoff is not None or indel_bin_cutoff is not None) and (aggregated_bin_ht is None or bin_id is None): 
        raise ValueError("If using bin cutoffs, aggregated_bin_ht and bin_id must be supplied.")

    # --- Step 4: Determine score cutoffs from bins if needed ---
    if snp_bin_cutoff is not None:
        snp_score_cutoff = aggregated_bin_ht.aggregate(
            hl.agg.filter(
                aggregated_bin_ht.snv & (aggregated_bin_ht.bin_id == bin_id) & (aggregated_bin_ht.bin == snp_bin_cutoff),
                hl.agg.min(aggregated_bin_ht.min_score),
            )
        )

    if indel_bin_cutoff is not None:
        indel_score_cutoff = aggregated_bin_ht.aggregate(
            hl.agg.filter(
                ~aggregated_bin_ht.snv & (aggregated_bin_ht.bin_id == bin_id) & (aggregated_bin_ht.bin == indel_bin_cutoff),
                hl.agg.min(aggregated_bin_ht.min_score),
            )
        )

    # A bin with no scored variants has no min_score; report that clearly.
    for kind, bin_cutoff, score_cutoff in (
        ("SNP", snp_bin_cutoff, snp_score_cutoff),
        ("indel", indel_bin_cutoff, indel_score_cutoff),
    ):
        if bin_cutoff is not None and score_cutoff is None:
            raise ValueError(
                f"{kind} bin {bin_cutoff} (bin_id '{bin_id}') contains no scored "
                f"variants, so no score cutoff can be derived from it. Check that "
                f"the RF scored {kind}s -- e.g. that annotate_for_rf was run with "
                f"--impute_features."
            )

    # --- Step 5: Determine global score range ---
    min_score = ht.aggregate(hl.agg.min(ht.score))
    max_score = ht.aggregate(hl.agg.max(ht.score))

    if snp_score_cutoff is not None:
        if not (min_score <= snp_score_cutoff <= max_score):
            raise ValueError("snp_score_cutoff is outside the score range of the table.")
        snp_cutoff_global = hl.struct(min_score=snp_score_cutoff)
    if indel_score_cutoff is not None:
        if not (min_score <= indel_score_cutoff <= max_score):
            raise ValueError("indel_score_cutoff is outside the score range of the table.")
        indel_cutoff_global = hl.struct(min_score=indel_score_cutoff)

    # --- Step 6: Build filter logic ---
    filters = {
        model_name: (
            hl.is_missing(ht.score)
            | (hl.is_snp(ht.alleles[0], ht.alleles[1]) & (ht.score < snp_cutoff_global.min_score))
            | (~hl.is_snp(ht.alleles[0], ht.alleles[1]) & (ht.score < indel_cutoff_global.min_score))
        ),
        "InbreedingCoeff": ht.inbreeding_coeff_flag,
        "AC0": ht.ac0_filter_expr,
    }

    # Ensure no missing values
    filters = {k: hl.or_else(v, False) for k, v in filters.items()}

    # --- Step 7: Annotation expressions ---
    annotations_expr = {}
    if model_name == "RF":
        annotations_expr.update(
            {
                "positive_train_site": hl.or_else(ht.positive_train_site, False),
                "rf_tp_probability": ht.rf_probability["TP"],
            }
        )

    if "feature_imputed" in ht.row:
        annotations_expr.update(
            {x: hl.or_missing(~ht.feature_imputed[x], ht[x]) for x in ht.row.feature_imputed}
        )

    annotations_expr.update(
        {
            "transmitted_singleton": hl.or_missing(ht.transmitted_singleton_from_ts, ht.transmitted_singleton)
        }
    )

    # --- Step 8: Final transmute ---
    ht = ht.transmute(
        filters=add_filters_expr(filters=filters),
        monoallelic=ht.monoallelic_flag_expr,
        ts_ac_filter=ht.ts_ac_filter_expr,
        **({score_name: ht.score} if score_name else {}),
        **annotations_expr,
    )

    # --- Step 9: Rename bin columns if present ---
    bin_names = [x for x in ht.row if x.endswith("bin")]
    bin_names = [
        (x, (x.split("adj_")[0] + x.split("adj_")[1] if len(x.split("adj_")) == 2 else "raw_" + x))
        for x in bin_names
    ]
    if bin_names:
        ht = ht.transmute(**{j: ht[i] for i, j in bin_names})

    # --- Step 10: Annotate globals ---
    if "bin_stats" in ht.row:
        bin_stats_struct = hl.struct(**{j: ht.bin_stats[i] for i, j in bin_names})
    else:
        bin_stats_struct = hl.struct()

    ht = ht.annotate_globals(
        bin_stats=bin_stats_struct,
        filtering_model=hl.struct(
            model_name=model_name,
            score_name=score_name,
            snv_cutoff=snp_cutoff_global,
            indel_cutoff=indel_cutoff_global,
        ),
        inbreeding_coeff_cutoff=inbreeding_coeff_cutoff,
    )

    # --- Step 11: VQSR annotation if provided ---
    if vqsr_ht:
        vqsr = vqsr_ht[ht.key]
        ht = ht.annotate(
            vqsr=hl.struct(
                AS_VQSLOD=vqsr.info.AS_VQSLOD,
                AS_culprit=vqsr.info.AS_culprit,
                NEGATIVE_TRAIN_SITE=vqsr.info.NEGATIVE_TRAIN_SITE,
                POSITIVE_TRAIN_SITE=vqsr.info.POSITIVE_TRAIN_SITE,
            ),
            SOR=vqsr.info.SOR,
        )

    #ht = ht.drop("AS_culprit") # Don't have this field
    return ht


def main(args):  # noqa: D103
    # Load the main score table
    ht = hl.read_table(os.path.join(BASE, WRITE_DIR, f"{args.model_id}_score_bins.ht"))

    # Filter centromere/telomere regions if requested
    if args.filter_centromere_telomere:
        cen_tel_ht = hl.read_table(
            os.path.join(BASE, "hail", "hg38.telomeresAndMergedCentromeres.ht")
        )
        ht = ht.filter(~hl.is_defined(cen_tel_ht[ht.locus]))

    # Filter low-quality sites
    info_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, "split_info.ht"))
    ht = ht.filter(~hl.or_else(info_ht[ht.key].AS_lowqual, False))

    # Drop 'info' column for VQSR models
    if args.model_id.startswith("vqsr_"):
        ht = ht.drop("info")

    # Load frequency HT for v3.1.2 subset release
    freq_path = os.path.join(BASE, WRITE_DIR, "corrected_freq.ht")
    release_path = os.path.join(BASE, WRITE_DIR, "release_sites.ht")
    if os.path.exists(freq_path):
        freq_ht = hl.read_table(freq_path)
    elif os.path.exists(release_path):
        freq_ht = hl.read_table(release_path)
    else:
        raise DataException(
            "There is no frequency HT or release sites HT available for the current release!"
        )

    # Select only the necessary fields from freq_ht
    freq_ht = freq_ht.select(
        InbreedingCoeff=freq_ht.InbreedingCoeff,
        freq_ac0=freq_ht.freq[0].AC,
        freq_monoallelic_af=freq_ht.freq[1].AF,
        freq_transmitted_singleton_ac=freq_ht.freq[1].AC,
    )

    # Join frequency info into ht
    freq_ht = freq_ht.key_by(*ht.key)
    ht = ht.join(freq_ht, how="left")

    # --- Annotate filters directly on ht (same table source) ---
    ht = ht.annotate(
        InbreedingCoeff_final=hl.coalesce(ht.InbreedingCoeff, ht.InbreedingCoeff_1),
        ac0_filter_expr=hl.or_else(ht.freq_ac0 == 0, False),
        monoallelic_flag_expr=hl.or_else((ht.freq_monoallelic_af == 1) | (ht.freq_monoallelic_af == 0), False),
        ts_ac_filter_expr=hl.or_else(ht.freq_transmitted_singleton_ac == 1, False)
    )

    # Load aggregated bin table
    aggregated_bin_path = os.path.join(BASE, WRITE_DIR, f"{args.model_id}_aggregated_bins.ht")
    if not file_exists(aggregated_bin_path):
        sys.exit(
            "Could not find binned HT for model:"
            f" {args.model_id} ({aggregated_bin_path}). Please run create_ranked_scores.py for that hash."
        )
    aggregated_bin_ht = hl.read_table(aggregated_bin_path)

    # Generate final filter HT using expressions from the same table
    ht = generate_final_filter_ht(
        ht,
        model_name=args.model_name,
        score_name=args.score_name,
        snp_bin_cutoff=args.snp_bin_cutoff,
        indel_bin_cutoff=args.indel_bin_cutoff,
        snp_score_cutoff=args.snp_score_cutoff,
        indel_score_cutoff=args.indel_score_cutoff,
        inbreeding_coeff_cutoff=args.inbreeding_coeff_threshold,
        aggregated_bin_ht=aggregated_bin_ht,
        bin_id="bin",
        vqsr_ht=None,
    )


    # Annotate globals
    ht = ht.annotate_globals(filtering_model=ht.filtering_model.annotate(model_id=args.model_id))

    if args.model_id.startswith("vqsr_"):
        ht = ht.annotate_globals(
            filtering_model=ht.filtering_model.annotate(
                snv_training_variables=[
                    "AS_QD", "AS_MQRankSum", "AS_ReadPosRankSum", "AS_FS", "AS_SOR", "AS_MQ",
                ],
                indel_training_variables=[
                    "AS_QD", "AS_MQRankSum", "AS_ReadPosRankSum", "AS_FS", "AS_SOR",
                ],
            )
        )
    else:
        ht = ht.annotate_globals(
            filtering_model=ht.filtering_model.annotate(
                snv_training_variables=ht.features,
                indel_training_variables=ht.features,
            )
        )

    ht.write(os.path.join(BASE, WRITE_DIR, "final_filter.ht"), args.overwrite)

    final_filter_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, "final_filter.ht"))

    n_filtered = hl.agg.count_where(
        (hl.len(final_filter_ht.filters) == 0) &
        (~final_filter_ht.fail_hard_filters)
    )

    summary = final_filter_ht.aggregate(hl.struct(
        n=hl.agg.count(),
        n_filtered=n_filtered,
        min_af=hl.agg.min(final_filter_ht.freq_monoallelic_af),
        max_af=hl.agg.max(final_filter_ht.freq_monoallelic_af),
    ))

    print(summary)
    



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset (default: False).",
        action="store_true",
    )
    parser.add_argument("--model_id", help="Filtering model ID to use.")
    parser.add_argument(
        "--model_name",
        help=(
            "Filtering model name to use in the filters field ('VQSR', 'AS_VQSR', or"
            " 'RF')."
        ),
        choices=["AS_VQSR", "VQSR", "RF"],
    )
    parser.add_argument(
        "--score_name",
        help=(
            "What to rename the filtering score annotation. This will be used in place"
            " of 'score' in the release HT info struct and the INFO field of the VCF"
            " (e.g. 'RF', 'AS_VQSLOD')."
        ),
    )
    parser.add_argument(
        "--inbreeding_coeff_threshold",
        help="InbreedingCoeff hard filter to use for variants.",
        type=float,
        default=INBREEDING_COEFF_HARD_CUTOFF,
    )
    snp_cutoff = parser.add_mutually_exclusive_group(required=True)
    snp_cutoff.add_argument(
        "--snp_bin_cutoff",
        help=(
            "RF or VQSR score bin to use as cutoff for SNPs. Value should be between 1"
            " and 100."
        ),
        type=float,
    )
    snp_cutoff.add_argument(
        "--snp_score_cutoff",
        help="RF or VQSR score to use as cutoff for SNPs.",
        type=float,
    )
    indel_cutoff = parser.add_mutually_exclusive_group(required=True)
    indel_cutoff.add_argument(
        "--indel_bin_cutoff",
        help=(
            "RF or VQSR score bin to use as cutoff for indels. Value should be between"
            " 1 and 100."
        ),
        type=float,
    )
    indel_cutoff.add_argument(
        "--indel_score_cutoff",
        help="RF or VQSR score to use as cutoff for indels.",
        type=float,
    )
    parser.add_argument(
        "--filter_centromere_telomere",
        help="Filter centromeres and telomeres from final filter Table.",
        action="store_true",
    )
    parser.add_argument(  # NOTE: This was required for v3.1 to grab the SOR annotation, we now compute this in `get_site_info_expr`
        "--vqsr_model_id",
        help=(
            "If a VQSR model ID is provided, a 'vqsr' annotation will be added to the"
            " final filter Table containing AS_VQSLOD , AS_culprit, NEGATIVE_TRAIN_SITE"
            " and POSITIVE_TRAIN_SITE."
        ),
        default="vqsr_alleleSpecificTrans",
        choices=["vqsr_classic", "vqsr_alleleSpecific", "vqsr_alleleSpecificTrans"],
        type=str,
    )
    parser.add_argument(
        "--hgdp_tgp_subset",
        help=(
            "Use HGDP + 1KG/TGP score binned filter Table instead of the full gnomAD"
            " Table."
        ),
        action="store_true",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
