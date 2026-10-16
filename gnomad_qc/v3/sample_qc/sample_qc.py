# %%
import pandas as pd
import argparse
import logging
import math
import pickle
from typing import Any, List, Tuple
import hail as hl
import os

from gnomad_qc.pipeline_config import BASE, WRITE_DIR, init_hail

init_hail("sample_qc")

from gnomad.assessment.validity_checks import compare_row_counts
from gnomad.resources.grch38.reference_data import (
    clinvar,
    lcr_intervals,
    purcell_5k_intervals,
    telomeres_and_centromeres,
)
from gnomad_qc.v3.resources.ancestry import assign_population_pcs, run_pca_with_relateds
from gnomad.sample_qc.filtering import (
    compute_qc_metrics_residuals,
    compute_stratified_metrics_filter,
    compute_stratified_sample_qc,
)
from gnomad_qc.v3.resources.pipeline import annotate_sex, get_qc_mt
from gnomad_qc.v3.resources.relatedness import (
    DUPLICATE_OR_TWINS,
    PARENT_CHILD,
    SIBLINGS,
    UNRELATED,
    compute_related_samples_to_drop,
    get_relationship_expr,
)
from gnomad.sample_qc.sex import get_ploidy_cutoffs, get_sex_expr
from gnomad.utils.annotations import bi_allelic_expr, get_adj_expr
from gnomad_qc.v3.resources.filtering import (
    add_filters_expr,
    filter_low_conf_regions,
    filter_to_autosomes,
    filter_to_clinvar_pathogenic,
)
from gnomad.utils.sparse_mt import densify_sites, filter_ref_blocks

from gnomad_qc.v2.resources.sample_qc import get_liftover_v2_qc_mt
from gnomad_qc.v3.resources.annotations import freq, get_info, last_END_position
from gnomad_qc.v3.resources.basics import get_gnomad_v3_mt
from gnomad_qc.v3.resources.meta import meta, meta_tsv_path, project_meta
from gnomad_qc.v3.resources.sample_qc import (
    ancestry_pca_eigenvalues,
    ancestry_pca_loadings,
    ancestry_pca_scores,
    get_sample_qc,
    hard_filtered_samples,
    pc_relate_pca_scores,
    pca_related_samples_to_drop,
    pca_samples_rankings,
    picard_metrics,
    pop,
    pop_rf_path,
    pop_tsv_path,
    qc,
    regressed_metrics,
    relatedness,
    release_related_samples_to_drop,
    release_samples_rankings,
    sample_clinvar_count,
    sample_inbreeding,
    sex,
    stratified_metrics,
)
# %%
logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("sample_qc")
logger.setLevel(logging.INFO)

SUBSETS = [
    "non_topmed",
    "controls_and_biobanks",
    "non_neuro",
    "non_v2",
    "non_cancer",
    "tgp",
    "hgdp",
]
"""
Default list of subsets to include in the sample QC meta HT from the project meta HT
"""


def compute_sample_qc() -> hl.Table:
    """
    Perform sample QC on the raw split matrix table using `compute_stratified_sample_qc`.

    :return: Table containing sample QC metrics
    :rtype: hl.Table
    """
    logger.info("Computing sample QC")

    
    mt = filter_to_autosomes(
        get_gnomad_v3_mt(
            split=True,
            key_by_locus_and_alleles=True,
            remove_hard_filtered_samples=False,
        )
    )
    mt = mt.select_entries(GT=mt.LGT)

    # Filter reference blocks
    mt = filter_ref_blocks(mt)

    # Remove centromeres and telomeres incase they were included
    mt = filter_low_conf_regions(
        mt,
        filter_lcr=False,
        filter_decoy=False,
        filter_segdup=False,
        filter_telomeres_and_centromeres=True,
    )

    sample_qc_ht = compute_stratified_sample_qc(
        mt,
        strata={
            "bi_allelic": bi_allelic_expr(mt),
            "multi_allelic": ~bi_allelic_expr(mt),
        },
        tmp_ht_prefix=os.path.join(BASE, WRITE_DIR, 'tmp_sample_qc.ht')
    )

    # Remove annotations that cannot be computed from the sparse format
    sample_qc_ht = sample_qc_ht.annotate(
        **{
            x: sample_qc_ht[x].drop(
                "n_called", "n_not_called", "n_filtered", "call_rate"
            )
            for x in sample_qc_ht.row_value
        }
    )

    return sample_qc_ht.repartition(100)


def compute_qc_mt(
    min_af: float = 0.0, min_inbreeding_coeff_threshold: float = -0.025
) -> hl.MatrixTable:
    """
    Generate the QC MT to be used in downstream sample QC.

    Filter MT to bi-allelic SNVs at gnomAD v2 QC and p5k sites. Compute and return QC MT.

    :param min_af: Minimum variant allele frequency to retain variant in qc matrix table
    :param min_inbreeding_coeff_threshold: Minimum site inbreeding coefficient to keep
    :return: MatrixTable filtered to variants for sample QC
    """
    # Load v2 and p5k sites for QC
    v2_qc_sites = hl.read_table(os.path.join(BASE, 'gnomad.genomes.r3.0.sites.ht'))
    v2_qc_sites = v2_qc_sites.key_by("locus")
    purcell5k = hl.read_table(os.path.join(BASE,'purcell5k.ht'))
    qc_sites = v2_qc_sites.union(purcell5k, unify=True)
    lcr_intervals = hl.read_table(os.path.join(BASE, 'hail', 'LCRFromHengHg38.ht'))
    qc_sites = v2_qc_sites.filter(hl.is_missing(lcr_intervals[v2_qc_sites.key]))

    mt = get_gnomad_v3_mt(
        key_by_locus_and_alleles=True,
        remove_hard_filtered_samples=False,
    )

    mt = mt.select_entries(
        "END", GT=mt.LGT, adj=get_adj_expr(mt.LGT, mt.GQ, mt.DP, mt.LAD)
    )
    mt = densify_sites(mt, qc_sites, hl.read_table(os.path.join(BASE, WRITE_DIR, 'HostSeq_last_END_positions.ht')))

    # Filter MT to bi-allelic SNVs that are found in the v2 QC and p5k HT
    mt = mt.filter_rows(
        (hl.len(mt.alleles) == 2)
        & hl.is_snp(mt.alleles[0], mt.alleles[1])
        & (qc_sites[mt.locus].alleles == mt.alleles)
    )
    mt = mt.checkpoint(
        os.path.join(BASE, WRITE_DIR, 'qc_mt_checkpoint.mt'), overwrite=True
    )
    mt = mt.naive_coalesce(5000)
    mt = mt.checkpoint(
        os.path.join(BASE, WRITE_DIR, 'qc_mt_checkpoint_repartitioned.mt'),
        overwrite=True,
    )

    info_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, 'split_info.ht'))
    info_ht = info_ht.annotate(
        info=info_ht.info.select(
            # No need for AS_annotations since it's bi-allelic sites only
            **{x: info_ht.info[x] for x in info_ht.info if not x.startswith("AS_")}
        )
    )
    mt = mt.annotate_rows(info=info_ht[mt.row_key].info)

    qc_mt = get_qc_mt(
        mt,
        min_af=min_af,
        min_inbreeding_coeff_threshold=min_inbreeding_coeff_threshold,
        min_hardy_weinberg_threshold=None,
        ld_r2=None,
        filter_lcr=False,
        filter_decoy=False,
        filter_segdup=False,
    )
    return qc_mt


def compute_hard_filters(
    cov_threshold: int = 10,
    max_n_snp: float = 2e6,
    min_n_snp: float = 5e5,
    max_n_singleton: float = 5e5,
    max_r_het_hom_var: float = 11,
    max_pct_contamination: float = 5.00,
    max_pct_chimera: float = 5.00,
    min_median_insert_size: int = 250,
    include_sex_filter: bool = True,
) -> hl.Table:
    """
    ht = hl.read_table(os.path.join(BASE, WRITE_DIR, 'tmp_sample_qc.ht_bi_allelic.ht'))
    summary = ht.aggregate(
        hl.struct(**{
            m: hl.struct(
                min=hl.agg.min(ht.sample_qc[m]),
                max=hl.agg.max(ht.sample_qc[m]),
                mean=hl.agg.mean(ht.sample_qc[m]),
                p5=hl.agg.approx_quantiles(ht.sample_qc[m], 0.05),
                p95=hl.agg.approx_quantiles(ht.sample_qc[m], 0.95),
            ) for m in metrics
        })
    )
    Struct(n_snp=Struct(min=738130, max=1792875, mean=969280.1596262038, p5=889552, p95=1501583), 
           n_singleton=Struct(min=0, max=76377, mean=9731.223800896349, p5=2630, p95=19414), 
           r_het_hom_var=Struct(min=1.8163662107050798, max=10.89423298351811, mean=3.4663717956456233, p5=2.8385344088022624, p95=4.910841708294667))

    Apply hard filters to samples and return Table with samples and the reason for filtering.

    This function expects a sex table generated by `impute_sex`.

    :param cov_threshold: Filtering threshold to use for chr20 coverage
    :param max_n_snp: Filtering threshold to use for the max number of SNPs
    :param min_n_snp: Filtering threshold to use for the min number of SNPs
    :param max_n_singleton: Filtering threshold to use for the max number of singletons
    :param max_r_het_hom_var: Filtering threshold to use for the max ratio of heterozygotes to alternate homozygotes
    :param max_pct_contamination: Filtering threshold to use for max percent contamination (this is a percent not a
        proportion, e.g. 5% == 5.00, %5 != 0.05)
    :param max_pct_chimera: Filtering threshold to use for max percent chimera (this is a percent not a proportion,
        e.g. 5% == 5.00, %5 != 0.05)
    :param min_median_insert_size: Filtering threshold to use for min median insert size
    :param include_sex_filter: Should sex inference be used in filtering
    :return: Table of hard filtered samples
    :rtype: hl.Table
    """
    ht = get_gnomad_v3_mt(remove_hard_filtered_samples=False).cols()
    hard_filters = dict()

    # Remove samples failing fingerprinting
    # TODO: Add these into hard filtering metadata when incorporating internal samples Picard metrics # noqa
    hard_filters["failed_fingerprinting"] = hl.array(
        ["09C90823", "10C103592", "S5530"]
    ).contains(ht.s)

    # Remove TCGA tumor samples based on TCGA naming convention:
    # https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
    hard_filters["TCGA_tumor_sample"] = ht.s.startswith("TCGA") & (
        hl.int(hl.str(ht.s).split("-")[3][:2]) < 10
    )

    # Remove low-coverage samples
    # chrom 20 coverage is computed to infer sex and used here
    cov_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, 'impute_sex.ht'))
    hard_filters["low_coverage"] = cov_ht[ht.key].chr20_mean_dp < cov_threshold

    # Remove extreme raw bi-allelic sample QC outliers
    # These were determined by visual inspection of the metrics
    
    bi_allelic_qc_struct = hl.read_table(os.path.join(BASE, WRITE_DIR, 'tmp_sample_qc.ht_bi_allelic.ht'))
    bi_allelic_qc_struct = bi_allelic_qc_struct[ht.key]
    hard_filters["bad_qc_metrics"] = (
        (bi_allelic_qc_struct.sample_qc.n_snp > max_n_snp)
        | (bi_allelic_qc_struct.sample_qc.n_snp < min_n_snp)
        | (bi_allelic_qc_struct.sample_qc.n_singleton > max_n_singleton)
        | (bi_allelic_qc_struct.sample_qc.r_het_hom_var > max_r_het_hom_var)
    )

    # Remove samples that fail picard metric thresholds
    ''''
    picard_ht = picard_metrics.ht()[ht.key]
    hard_filters["contamination"] = (
        picard_ht.bam_metrics.freemix > max_pct_contamination
    )
    hard_filters["chimera"] = picard_ht.bam_metrics.pct_chimeras > max_pct_chimera
    hard_filters["insert_size"] = (
        picard_ht.bam_metrics.median_insert_size < min_median_insert_size
    )
    '''

    if include_sex_filter:
        # Remove samples with ambiguous sex assignments
        sex_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, 'impute_sex.ht'))
        sex_ht = sex_ht[ht.key]
        hard_filters["ambiguous_sex"] = sex_ht.sex_karyotype == "ambiguous"
        hard_filters["sex_aneuploidy"] = ~hl.set({"ambiguous", "XX", "XY"}).contains(
            sex_ht.sex_karyotype
        )

    ht = ht.annotate(hard_filters=add_filters_expr(filters=hard_filters))

    ht = ht.filter(hl.len(ht.hard_filters) > 0)
    ht = ht.annotate_globals(
        hard_filter_cutoffs=hl.struct(
            min_cov=cov_threshold,
            max_n_snp=max_n_snp,
            min_n_snp=min_n_snp,
            max_n_singleton=max_n_singleton,
            max_r_het_hom_var=max_r_het_hom_var,
            max_pct_contamination=max_pct_contamination,
            max_pct_chimera=max_pct_chimera,
            min_median_insert_size=min_median_insert_size,
        ),
    )

    # Drop extra MT columns that shouldn't be in the hard filter table
    cols_to_drop = [c for c in ht.row if c not in ['s', 'hard_filters']]
    if cols_to_drop:
        ht = ht.drop(*cols_to_drop)

    return ht


def compute_sex(aaf_threshold: float = 0.001, f_stat_cutoff: float = 0.5) -> hl.Table:
    """
    Impute sample sex based on X-chromosome heterozygosity and sex chromosome ploidy.

    Allele frequencies from v3.0 are used to prevent the need to densify the current version's sparse MT.

    :param aaf_threshold: Minimum alternate allele frequency to be used in f-stat calculations.
    :param f_stat_cutoff: f-stat to roughly divide 'XX' from 'XY' samples. Assumes XX samples are below cutoff and XY are above cutoff.
    :return: Table with inferred sex annotation
    :rtype: hl.Table
    """
    mt_pass = get_gnomad_v3_mt(
        key_by_locus_and_alleles=True,
        remove_hard_filtered_samples=False,
    )
    mt_pass = mt_pass.annotate_entries(GT=hl.vds.lgt_to_gt(mt_pass.LGT, mt_pass.LA))

    # Use AF from v3
    mt = hl.experimental.sparse_split_multi(mt_pass, filter_changed_loci=True)
    mt = hl.experimental.densify(mt)
    mt = hl.variant_qc(mt)

    freq_ht = mt.rows().select("variant_qc")
    freq_ht = freq_ht.filter(hl.len(freq_ht.variant_qc.AF) > 1)
    freq_ht = freq_ht.select(AF=freq_ht.variant_qc.AF[1])
    freq_ht = freq_ht.filter(freq_ht.AF > aaf_threshold)

    sex_ht = annotate_sex(
        mt_pass,
        excluded_intervals=hl.read_table(os.path.join(BASE, 'hail', 'hg38.telomeresAndMergedCentromeres.ht')),
        aaf_threshold=aaf_threshold,
        f_stat_cutoff=f_stat_cutoff,
        sites_ht=freq_ht,
        aaf_expr="AF",
        gt_expr="GT",
        use_gaussian_mixture_model=True,
    )

    return sex_ht


def reannotate_sex(
    cov_threshold: int,
    x_ploidy_cutoffs: Tuple[float, Tuple[float, float], float],
    y_ploidy_cutoffs: Tuple[Tuple[float, float], float],
):
    """
    Rerun sex karyotyping annotations without re-computing sex imputation metrics.

    :param cov_threshold: Filtering threshold to use for chr20 coverage in `compute_hard_filters`
    :param x_ploidy_cutoffs: Tuple of X chromosome ploidy cutoffs: (upper cutoff for single X, (lower cutoff for double X, upper cutoff for double X), lower cutoff for triple X)
    :param y_ploidy_cutoffs: Tuple of Y chromosome ploidy cutoffs: ((lower cutoff for single Y, upper cutoff for single Y), lower cutoff for double Y)
    :return: Table with sex karyotyping annotations
    :rtype: hl.Table
    """
    # Copy HT to temp location to overwrite annotation
    sex_ht = sex.ht().checkpoint("gs://gnomad-tmp/sex_ht_checkpoint.ht", overwrite=True)
    hard_filter_ht = compute_hard_filters(cov_threshold, include_sex_filter=False)

    # Copy HT to temp location because it uses sex_ht for chr20 coverage
    hard_filter_ht = hard_filter_ht.checkpoint(
        "gs://gnomad-tmp/hardfilter_checkpoint.ht", overwrite=True
    )
    new_x_ploidy_cutoffs, new_y_ploidy_cutoffs = get_ploidy_cutoffs(
        sex_ht.filter(hl.is_missing(hard_filter_ht[sex_ht.key])), f_stat_cutoff=0.5
    )
    x_ploidy_cutoffs = hl.struct(
        upper_x=x_ploidy_cutoffs[0] if x_ploidy_cutoffs[0] else new_x_ploidy_cutoffs[0],
        lower_xx=(
            x_ploidy_cutoffs[1][0]
            if x_ploidy_cutoffs[1][0]
            else new_x_ploidy_cutoffs[1][0]
        ),
        upper_xx=(
            x_ploidy_cutoffs[1][1]
            if x_ploidy_cutoffs[1][1]
            else new_x_ploidy_cutoffs[1][1]
        ),
        lower_xxx=(
            x_ploidy_cutoffs[2] if x_ploidy_cutoffs[2] else new_x_ploidy_cutoffs[2]
        ),
    )
    y_ploidy_cutoffs = hl.struct(
        lower_y=(
            y_ploidy_cutoffs[0][0]
            if y_ploidy_cutoffs[0][0]
            else new_y_ploidy_cutoffs[0][0]
        ),
        upper_y=(
            y_ploidy_cutoffs[0][1]
            if y_ploidy_cutoffs[0][1]
            else new_y_ploidy_cutoffs[0][1]
        ),
        lower_yy=(
            y_ploidy_cutoffs[1] if y_ploidy_cutoffs[1] else new_y_ploidy_cutoffs[1]
        ),
    )
    sex_ht = sex_ht.annotate(
        **get_sex_expr(
            sex_ht.chrX_ploidy,
            sex_ht.chrY_ploidy,
            (
                x_ploidy_cutoffs["upper_x"],
                (x_ploidy_cutoffs["lower_xx"], x_ploidy_cutoffs["upper_xx"]),
                x_ploidy_cutoffs["lower_xxx"],
            ),
            (
                (y_ploidy_cutoffs["lower_y"], y_ploidy_cutoffs["upper_y"]),
                y_ploidy_cutoffs["lower_yy"],
            ),
        )
    )
    sex_ht = sex_ht.annotate_globals(
        x_ploidy_cutoffs=x_ploidy_cutoffs,
        y_ploidy_cutoffs=y_ploidy_cutoffs,
    )

    return sex_ht


def compute_sample_rankings(use_qc_metrics_filters: bool) -> hl.Table:
    """
    Rank samples by sample filtering and chr20 mean coverage.

    Filtered samples are ranked lowest, followed by non-releasable samples, and then by coverage where higher coverage
    is a higher rank.

    :param use_qc_metrics_filters: Should sample QC metric filters be included or only hard filters
    :return: Table with sample rankings
    :rtype: hl.Table
    """
    sex_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, 'impute_sex.ht')).key_by("s")
    project_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, 'meta.ht')).key_by("s")

    project_ht = project_ht.select(
        "releasable",
        "exclude",
        chr20_mean_dp=sex_ht[project_ht.key].chr20_mean_dp,
        filtered=False, # No hard filtered samples
    )

    if use_qc_metrics_filters:
        project_ht = project_ht.annotate(
            filtered=hl.cond(
                project_ht.filtered,
                True,
                False,
            )
        )
    project_ht.describe()
    project_ht = project_ht.order_by(
        project_ht.filtered,
        hl.desc(project_ht.releasable & ~project_ht.exclude),
        hl.desc(project_ht.chr20_mean_dp),
    ).add_index(name="rank")
    project_ht.describe()
    return project_ht.key_by("s").select("filtered", "rank")


def run_pca(
    include_unreleasable_samples: bool, n_pcs: int, related_samples_to_drop: hl.Table
) -> Tuple[List[float], hl.Table, hl.Table]:
    """
    Run population PCA using `run_pca_with_relateds`.

    :param include_unreleasable_samples: Should unreleasable samples be included in the PCA
    :param n_pcs: Number of PCs to compute
    :param related_samples_to_drop: Table of related samples to drop from PCA run
    :return: eigenvalues, scores and loadings
    """
    logger.info("Running population PCA")
    qc_mt = hl.read_matrix_table(os.path.join(BASE, WRITE_DIR, 'compute_qc_mt.mt'))
    # Read sub-sampled table from disk
    sub_sample_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, 'sub_sampled_qc.ht'))
    # Use semi_join_rows to get filtered mt to feed to hwe_normalized_pca
    filtered_qc_mt = qc_mt.semi_join_rows(sub_sample_ht)

    #samples_to_drop = related_samples_to_drop.select()
    samples_to_drop = hl.read_table(os.path.join(BASE, WRITE_DIR, 'samples_to_drop.ht'))
    if not include_unreleasable_samples:
        logger.info("Excluding unreleasable samples for PCA.")
        # TODO: Fix keys/re-enable this for future use (don't need it for the time being)
        '''
        samples_to_drop = samples_to_drop.union(
            qc_mt.filter_cols(
                ~project_meta.ht()[qc_mt.col_key].releasable
                | project_meta.ht()[qc_mt.col_key].exclude
            )
            .cols()
            .select()
        )
        '''
    else:
        logger.info("Including unreleasable samples for PCA")

    return run_pca_with_relateds(filtered_qc_mt, samples_to_drop, n_pcs=n_pcs)


def calculate_mislabeled_training(pop_ht: hl.Table, pop_field: str):# -> [int, float]:
    """
    Calculate the number and proportion of mislabeled training samples.

    :param pop_ht: Table with assigned pops/subpops that is returned by `assign_population_pcs`
    :param pop_field: Name of field in the Table containing the assigned pop/subpop
    :return: The number and proportion of mislabeled training samples
    """
    n_mislabeled_samples = pop_ht.aggregate(
        hl.agg.count_where(pop_ht.training_pop != pop_ht[pop_field])
    )

    defined_training_pops = pop_ht.aggregate(
        hl.agg.count_where(hl.is_defined(pop_ht.training_pop))
    )

    prop_mislabeled_samples = n_mislabeled_samples / defined_training_pops

    return n_mislabeled_samples, prop_mislabeled_samples


def assign_pops(
    min_prob: float,
    include_unreleasable_samples: bool,
    max_number_mislabeled_training_samples: int = None,
    max_proportion_mislabeled_training_samples: float = None,
    pcs: List[int] = list(range(1, 17)),
    withhold_prop: float = None,
    pop: str = None,
    curated_subpops: list = None,
    additional_samples_to_drop: hl.Table = None,
    high_quality: bool = False,
    missing_label: str = "oth",
    seed: int = 24,
) -> Tuple[hl.Table, Any]:
    """
    Use a random forest model to assign global population labels based on the results from `assign_pca`.

    The training data used is `project_pop` on the project metadata HT where it is defined, otherwise it uses the
    `v2_pop` if defined and not `oth`.

    The method starts by inferring the pop on all samples, then comparing the training data to the inferred data,
    removing the truth outliers and re-training. This is repeated until the number of truth samples is less than
    or equal to `max_number_mislabeled_training_samples` or `max_proportion_mislabeled_training_samples`. Only one
    of `max_number_mislabeled_training_samples` or `max_proportion_mislabeled_training_samples` can be set.

    :param min_prob: Minimum RF probability for pop assignment
    :param include_unreleasable_samples: Should unreleasable samples be included in the PCA
    :param max_number_mislabeled_training_samples: Maximum number of training samples that can be mislabelled
    :param max_proportion_mislabeled_training_samples: Maximum proportion of training samples that can be mislabelled
    :param pcs: List of PCs to use in the RF
    :param withhold_prop: Proportion of training pop samples to withhold from training will keep all samples if `None`
    :param pop: Population that the PCA was restricted to. When set to None, the PCA on the full dataset is returned
    :param curated_subpops: List of the subpops that are accepted as known labels for the pop. If not specified, all values of subpop_description will be used as known labels when inferring subpops.
    :param additional_samples_to_drop: Optional table of additional samples to drop when training the random forest model
    :param high_quality: Whether the file includes PCA info for only high-quality samples
    :param missing_label: Label for samples for which the assignment probability is smaller than `min_prob`
    :param seed: Random seed
    :return: Table of pop or subpop assignments and the RF model
    :rtype: hl.Table
    """
    logger.info("Assigning global population/subpopulation labels")
    pop_pca_scores_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, 'ancestry_pca_scores.ht'))
    project_meta_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, "meta.ht"))
    #project_meta_ht = project_meta[pop_pca_scores_ht["s"]]
    
    if additional_samples_to_drop is not None:
        samples_to_drop = hl.literal(additional_samples_to_drop.s.collect())

    if pop is not None:
        # Set subpop_description to missing if the subpop is not in the manually
        # curated list of subpops for the given pop
        if curated_subpops is not None:
            pop_pca_scores_ht = pop_pca_scores_ht.annotate(
                subpop_description=hl.or_missing(
                    hl.literal(curated_subpops).contains(
                        project_meta_ht.subpop_description
                    ),
                    project_meta_ht.subpop_description,
                )
            )
        else:
            pop_pca_scores_ht = pop_pca_scores_ht.annotate(
                subpop_description=project_meta_ht.subpop_description
            )

        if additional_samples_to_drop is not None:
            pop_pca_scores_ht = pop_pca_scores_ht.annotate(
                subpop_description=hl.or_missing(
                    ~samples_to_drop.contains(pop_pca_scores_ht.s),
                    pop_pca_scores_ht.subpop_description,
                )
            )

        pop_pca_scores_ht = pop_pca_scores_ht.annotate(
            training_pop=(
                hl.case()
                .when(
                    hl.is_defined(pop_pca_scores_ht.subpop_description),
                    pop_pca_scores_ht.subpop_description,
                )
                .or_missing()
            )
        )

        pop_field = "subpop"
    else:
        pop_pca_scores_ht = pop_pca_scores_ht.annotate(
            training_pop=True
        )

        if additional_samples_to_drop is not None:
            pop_pca_scores_ht = pop_pca_scores_ht.annotate(
                training_pop=hl.or_missing(
                    ~samples_to_drop.contains(pop_pca_scores_ht.s),
                    pop_pca_scores_ht.training_pop,
                )
            )

        # Keep track of defined pops, useful if samples are withheld
        pop_pca_scores_ht = pop_pca_scores_ht.annotate(
            defined_pop=pop_pca_scores_ht.training_pop
        )
        pop_field = "pop"

    pop_pca_scores_ht = pop_pca_scores_ht.annotate(
        training_pop_all=pop_pca_scores_ht.training_pop
    )

    if withhold_prop:
        pop_pca_scores_ht = pop_pca_scores_ht.annotate(
            training_pop=hl.or_missing(
                hl.is_defined(pop_pca_scores_ht.training_pop)
                & hl.rand_bool(1.0 - withhold_prop, seed=seed),
                pop_pca_scores_ht.training_pop,
            )
        )

    if pop is not None:
        pop_pca_scores_ht = pop_pca_scores_ht.annotate(
            withheld_sample=hl.is_defined(pop_pca_scores_ht.subpop_description)
            & (~hl.is_defined(pop_pca_scores_ht.training_pop))
        )
    else:
        pop_pca_scores_ht = pop_pca_scores_ht.annotate(
            withheld_sample=hl.is_defined(pop_pca_scores_ht.defined_pop)
            & (~hl.is_defined(pop_pca_scores_ht.training_pop))
        )

    if (
        max_number_mislabeled_training_samples is not None
        and max_proportion_mislabeled_training_samples is None
    ):
        max_mislabeled = max_number_mislabeled_training_samples
    elif (
        max_proportion_mislabeled_training_samples is not None
        and max_number_mislabeled_training_samples is None
    ):
        max_mislabeled = max_proportion_mislabeled_training_samples
    else:
        raise ValueError(
            "One and only one of max_number_mislabeled_training_samples or"
            " max_proportion_mislabeled_training_samples must be set!"
        )

    logger.info(
        "Running RF using {} training examples".format(
            pop_pca_scores_ht.aggregate(
                hl.agg.count_where(hl.is_defined(pop_pca_scores_ht.training_pop))
            )
        )
    )

    pop_ht, pops_rf_model = assign_population_pcs(
        pop_pca_scores_ht,
        pc_cols=pcs,
        known_col="training_pop",
        output_col=pop_field,
        min_prob=min_prob,
        missing_label=missing_label,
    )

    # Calculate number and proportion of mislabeled samples
    n_mislabeled_samples, prop_mislabeled_samples = calculate_mislabeled_training(
        pop_ht, pop_field
    )

    if max_number_mislabeled_training_samples:
        mislabeled = n_mislabeled_samples
    else:
        mislabeled = prop_mislabeled_samples

    pop_ht = pop_ht.annotate(
        training_pop_all=pop_pca_scores_ht[pop_ht.key].training_pop_all
    )
    pop_assignment_iter = 1
    while mislabeled > max_mislabeled:
        pop_assignment_iter += 1
        logger.info(
            f"Found {n_mislabeled_samples}({round(prop_mislabeled_samples*100, 2)}%)"
            " training samples labeled differently from their known pop. Re-running"
            " without them."
        )

        pop_ht = pop_ht[pop_pca_scores_ht.key]
        pop_pca_scores_ht = pop_pca_scores_ht.annotate(
            training_pop=hl.or_missing(
                (pop_ht.training_pop == pop_ht[pop_field]),
                pop_pca_scores_ht.training_pop,
            ),
            training_pop_all=hl.or_missing(
                hl.is_missing(pop_ht.training_pop)
                | (pop_ht.training_pop == pop_ht[pop_field]),
                pop_pca_scores_ht.training_pop_all,
            ),
        ).persist()

        logger.info(
            "Running RF using {} training examples".format(
                pop_pca_scores_ht.aggregate(
                    hl.agg.count_where(hl.is_defined(pop_pca_scores_ht.training_pop))
                )
            )
        )

        pop_ht, pops_rf_model = assign_population_pcs(
            pop_pca_scores_ht,
            pc_cols=pcs,
            known_col="training_pop",
            output_col=pop_field,
            min_prob=min_prob,
            missing_label=missing_label,
        )

        pop_ht = pop_ht.annotate(
            training_pop_all=pop_pca_scores_ht[pop_ht.key].training_pop_all
        )

        # Calculate number and proportion of mislabeled samples
        n_mislabeled_samples, prop_mislabeled_samples = calculate_mislabeled_training(
            pop_ht, pop_field
        )

        if max_number_mislabeled_training_samples:
            mislabeled = n_mislabeled_samples
        else:
            mislabeled = prop_mislabeled_samples

    pop_ht = pop_ht.annotate_globals(
        min_prob=min_prob,
        include_unreleasable_samples=include_unreleasable_samples,
        max_mislabeled=max_mislabeled,
        pop_assignment_iterations=pop_assignment_iter,
        pcs=pcs,
        high_quality=high_quality,
        n_mislabeled_training_samples=n_mislabeled_samples,
        prop_mislabeled_training_samples=prop_mislabeled_samples,
    )
    # Only set when assigning subpopulations within `pop`.
    if pop is not None:
        pop_ht = pop_ht.annotate_globals(pop=pop)
    if withhold_prop:
        pop_ht = pop_ht.annotate_globals(withhold_prop=withhold_prop)
    pop_ht = pop_ht.annotate(
        withheld_sample=pop_pca_scores_ht[pop_ht.key].withheld_sample
    )

    return pop_ht, pops_rf_model


def apply_stratified_filters(
    sample_qc_ht: hl.Table, filtering_qc_metrics: List[str]
) -> hl.Table:
    """
    Use population stratified QC metrics to determine what samples are outliers and should be filtered.

    Use `compute_stratified_metrics_filter` to run hl.sample_qc and return the metrics stratified by assigned
    population with a `qc_metrics_filters` annotation indicating if the sample falls a certain number of MAD
    outside the distribution.

    .. note::

        Default left and right MAD for `compute_stratified_metrics_filter` is 4, we modify this for n_singleton to be
        (4.0, 8.0)

    :param sample_qc_ht: Sample QC HT
    :param filtering_qc_metrics: Specific metrics to compute
    :return: Table with stratified metrics and filters
    :rtype: hl.Table
    """
    logger.info(
        "Computing stratified QC metrics filters using metrics: "
        + ", ".join(filtering_qc_metrics)
    )
    pop_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, "ancestry_rf.assignment.legacy.ht"))
    hard_filtered_samples_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, 'compute_hard_filters.ht'))
    sample_qc_ht = sample_qc_ht.annotate(qc_pop=pop_ht[sample_qc_ht.key].pop)
    sample_qc_ht = sample_qc_ht.filter(
        hl.is_missing(hard_filtered_samples_ht[sample_qc_ht.key])
    )
    stratified_metrics_ht = compute_stratified_metrics_filter(
        sample_qc_ht,
        qc_metrics={
            metric: sample_qc_ht.sample_qc[metric] for metric in filtering_qc_metrics
        },
        strata={"qc_pop": sample_qc_ht.qc_pop},
        metric_threshold={"n_singleton": (4.0, 8.0)},
    )
    return stratified_metrics_ht


def apply_regressed_filters(
    sample_qc_ht: hl.Table,
    filtering_qc_metrics: List[str],
    include_unreleasable_samples: bool,
    n_pcs: int = 16,
) -> hl.Table:
    """
    Compute sample QC metrics residuals after regressing out population PCs and determine what samples are outliers that should be filtered.

    .. note::

        Default left and right MAD for `compute_stratified_metrics_filter` is 4. We modify this for n_singleton_residual
        to be (math.inf, 8.0) and for r_het_hom_var_residual to be (math.inf, 4.0). The math.inf is used to prevent a
        lower cutoff for these metrics.

    :param sample_qc_ht: Sample QC HT
    :param filtering_qc_metrics: Specific metrics to compute
    :param include_unreleasable_samples: Should unreleasable samples be included in the regression
    :param n_pcs: Number of population PCs to use in regression
    :return: Table with stratified metrics and filters
    :rtype: hl.Table
    """
    project_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, "meta.ht")) # project_meta.ht()
    project_ht = project_ht.key_by('s')
    ancestry_pca_scores = hl.read_table(os.path.join(BASE, WRITE_DIR, 'ancestry_pca_scores.ht'))
    hard_filtered_samples = hl.read_table(os.path.join(BASE, WRITE_DIR, 'compute_hard_filters.ht'))
    sample_qc_ht = sample_qc_ht.select(
        **sample_qc_ht.sample_qc,
        **ancestry_pca_scores[sample_qc_ht.key],
        releasable=project_ht[sample_qc_ht.key].releasable,
        exclude=project_ht[sample_qc_ht.key].exclude,
    )
    residuals_ht = compute_qc_metrics_residuals(
        ht=sample_qc_ht,
        pc_scores=sample_qc_ht.scores[:n_pcs],
        qc_metrics={metric: sample_qc_ht[metric] for metric in filtering_qc_metrics},
        regression_sample_inclusion_expr=sample_qc_ht.releasable
        & ~sample_qc_ht.exclude,
    )
    residuals_ht = residuals_ht.filter(
        hl.is_missing(hard_filtered_samples[residuals_ht.key])
    )
    stratified_metrics_ht = compute_stratified_metrics_filter(
        ht=residuals_ht,
        qc_metrics=dict(residuals_ht.row_value),
        metric_threshold={
            "n_singleton_residual": (math.inf, 8.0),
            "r_het_hom_var_residual": (math.inf, 4.0),
        },
    )

    residuals_ht = residuals_ht.annotate(**stratified_metrics_ht[residuals_ht.key])
    residuals_ht = residuals_ht.annotate_globals(
        **stratified_metrics_ht.index_globals(),
        n_pcs=n_pcs,
    )

    return residuals_ht


def get_relatedness_set_ht(relatedness_ht: hl.Table) -> hl.Table:
    """
    Parse relatedness Table to get every relationship (except UNRELATED) per sample.

    Return Table keyed by sample with all sample relationships in a set.

    :param Table relatedness_ht: Table with inferred relationship information output by pc_relate.
        Keyed by sample pair (i, j).
    :return: Table keyed by sample (s) with all relationships annotated as a set.
    :rtype: hl.Table
    """
    relatedness_ht = relatedness_ht.filter(relatedness_ht.relationship != UNRELATED)
    relatedness_ht = relatedness_ht.select("relationship", s=relatedness_ht.i.s).union(
        relatedness_ht.select("relationship", s=relatedness_ht.j.s)
    )
    relatedness_ht = relatedness_ht.group_by(relatedness_ht.s).aggregate(
        relationships=hl.agg.collect_as_set(relatedness_ht.relationship)
    )
    return relatedness_ht


def get_relationship_filter_expr(
    hard_filtered_expr: hl.expr.BooleanExpression,
    relationship: str,
    relationship_set: hl.expr.SetExpression,
) -> hl.expr.builders.CaseBuilder:
    """
    Return case statement to populate relatedness filters in sample_filters struct.

    :param hl.expr.BooleanExpression hard_filtered_expr: Boolean for whether sample was hard filtered.
    :param str relationship: Relationship to check for. One of DUPLICATE_OR_TWINS, PARENT_CHILD, or SIBLINGS.
    :param hl.expr.SetExpression relationship_set: Set containing all possible relationship strings for sample.
    :return: Case statement used to population sample_filters related filter field.
    :rtype: hl.expr.builders.CaseBuilder
    """
    return (
        hl.case()
        .when(hard_filtered_expr, hl.null(hl.tbool))
        .when(hl.is_defined(relationship_set), relationship_set.contains(relationship))
        .default(False)
    )


def join_tables(
    left_ht: hl.Table,
    left_key: str,
    right_ht: hl.Table,
    right_key: str,
    join_type: str,
    sample_count_match: bool = True,
) -> hl.Table:
    """
    Join left and right tables using specified keys and join types and returns result.

    Also prints warning if sample counts are not the same.

    :param Table left_ht: Left Table to be joined
    :param str left_key: Key of left Table
    :param Table right_ht: Right Table to be joined
    :param str right_key: Key of right Table
    :param str join_type: Type of join
    :param bool sample_count_match: Are the sample counts expected to match in the tables
    :return: Table with annotations
    :rtype: Table
    """
    if sample_count_match and not compare_row_counts(left_ht, right_ht):
        logger.warning("Sample counts in left and right tables do not match!")

        in_left_not_right = left_ht.anti_join(right_ht)
        if in_left_not_right.count() != 0:
            logger.warning(
                f"The following {in_left_not_right.count()} samples are found in the"
                " left HT, but are not found in the right HT"
            )
            in_left_not_right.select().show(n=-1)

        in_right_not_left = right_ht.anti_join(left_ht)
        if in_right_not_left.count() != 0:
            logger.warning(
                f"The following {in_right_not_left.count()} samples are found in the"
                " right HT, but are not found in left HT"
            )
            in_right_not_left.select().show(n=-1)

        if join_type != "outer":
            logger.warning(
                "Join type is not an outer join so some samples will be filtered!"
            )

    return left_ht.key_by(left_key).join(right_ht.key_by(right_key), how=join_type)


def generate_metadata(regressed_metrics_outlier: bool = True) -> hl.Table:
    """
    Pull all sample QC information together in one Table.

    :param regressed_metrics_outlier: Should the outlier table be from the residuals instead of pop stratified
    :return: Table annotated with all metadata
    :rtype: hl.Table
    """
    logging_statement = "Reading in {} and joining with meta HT"

    logger.info(
        "Loading metadata file with subset, age, and releasable information to begin"
        " creation of the meta HT"
    )
    left_ht = get_gnomad_v3_mt(remove_hard_filtered_samples=False).cols()

    right_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, "meta.ht")) # project_meta.ht()
    if '' in right_ht.row:
        right_ht = right_ht.drop('')
  
    right_ht = right_ht.select(
        s=right_ht.s,
        project_meta=hl.struct(**right_ht.row.drop("s")),
    ).key_by("s")

    left_ht = join_tables(left_ht, "s", right_ht.select_globals(), "s", "right")

    #logger.info(logging_statement.format("picard metric HT"))
    #right_ht = picard_metrics.ht()
    #right_ht = right_ht.select("bam_metrics")
    #left_ht = join_tables(left_ht, "s", right_ht, "s", "left", sample_count_match=False)

    logger.info(logging_statement.format("sex HT"))
    impute_stats = ["f_stat", "n_called", "expected_homs", "observed_homs"]
    
    right_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, "impute_sex.ht"))
    right_ht = right_ht.transmute(
        impute_sex_stats=hl.struct(**{x: right_ht[x] for x in impute_stats})
    )

    # Strip ref_block_max_length from globals before wrapping
    sex_globals_fields = [f for f in right_ht.globals.dtype.fields if f != 'ref_block_max_length']
    right_ht = right_ht.select_globals(
        sex_imputation_ploidy_cutoffs=hl.struct(**{f: right_ht.globals[f] for f in sex_globals_fields})
    )
    right_ht = right_ht.select(sex_imputation=hl.struct(**right_ht.row.drop("s")))
    left_ht = join_tables(left_ht, "s", right_ht, "s", "right")

    logger.info(logging_statement.format("sample QC HT"))
    right_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, "tmp_sample_qc.ht_bi_allelic.ht")) # get_sample_qc("bi_allelic").ht()

    # Drop extra MT fields before joining
    extra_fields = [f for f in right_ht.row if f not in ['s', 'sample_qc', 'bi_allelic', 'multi_allelic']]
    if extra_fields:
        right_ht = right_ht.drop(*extra_fields)
        
    right_ht = right_ht.select_globals()

    not_in_sparse = ["n_called", "n_not_called", "n_filtered", "call_rate"]
    right_ht = right_ht.annotate(
        **{
            x: right_ht[x].drop(*[f for f in not_in_sparse if f in right_ht[x].dtype.fields])
            for x in right_ht.row_value
            if isinstance(right_ht[x].dtype, hl.tstruct)
        }
    )

    left_ht = join_tables(left_ht, "s", right_ht, "s", "right")

    ancestry_bgz_ht = hl.import_table(
        os.path.join(BASE, WRITE_DIR, "ancestry_rf.assignment.tsv.bgz"),
        impute=True,
        types={"s": hl.tstr}
    )
    ancestry_bgz_ht = ancestry_bgz_ht.annotate(
        pca_scores=hl.parse_json(
            ancestry_bgz_ht.pca_scores,
            hl.tarray(hl.tfloat64)
        )
    )
    ancestry_bgz_ht = ancestry_bgz_ht.key_by("s")
    ancestry_bgz_ht.write(
        os.path.join(BASE, WRITE_DIR, "ancestry_rf.assignment.legacy.ht"),
        overwrite=True
    )

    logger.info(logging_statement.format("population PCA HT"))
    right_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, "ancestry_rf.assignment.legacy.ht")) # pop.ht()
    right_ht = right_ht.select_globals(
        population_inference_pca_metrics=right_ht.globals
    )
    right_ht = right_ht.select(population_inference=hl.struct(**right_ht.row.drop("s")))
    left_ht = join_tables(
        left_ht, "s", right_ht, "s", "outer", sample_count_match=False
    )

    logger.info(
        "Reading hard filters HT, renaming hard filters struct to sample_filters, and"
        " joining with meta HT"
    )
    right_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, "compute_hard_filters.ht")) # hard_filtered_samples.ht()

    # Drop extra row fields
    extra_fields = [f for f in right_ht.row if f not in ['s', 'hard_filters']]
    if extra_fields:
        right_ht = right_ht.drop(*extra_fields)

    # Strip ref_block_max_length from globals — only keep hard_filter_cutoffs
    right_ht = right_ht.select_globals('hard_filter_cutoffs')

    left_ht = join_tables(
        left_ht, "s", right_ht, "s", "outer", sample_count_match=False
    )
    left_ht = left_ht.checkpoint(os.path.join(BASE, WRITE_DIR, "population_joined.ht"), overwrite=True)

    # Change sample_filters to a struct
    ex_right_ht = right_ht.explode(right_ht.hard_filters)
    hard_filters = ex_right_ht.aggregate(
        hl.agg.collect_as_set(ex_right_ht.hard_filters)
    )
    left_ht = left_ht.transmute(
        sample_filters=hl.struct(
            **{
                v: hl.if_else(
                    hl.is_defined(left_ht.hard_filters),
                    left_ht.hard_filters.contains(v),
                    False,
                )
                for v in hard_filters
            },
            hard_filters=left_ht.hard_filters,
            hard_filtered=hl.is_defined(left_ht.hard_filters)
            & (hl.len(left_ht.hard_filters) > 0),
        )
    )
    left_ht = left_ht.checkpoint(os.path.join(BASE, WRITE_DIR, "sample_filters.ht"), overwrite=True)

    logger.info(
        "Reading in PCA related samples to drop HT and preparing to annotate meta HT's"
        " sample_filter struct with relatedness booleans"
    )
    related_samples_to_drop_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, 'samples_to_drop.ht')) # pca_related_samples_to_drop.ht()
    release_related_samples_to_drop_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, "release_related_samples_to_drop.ht")) # release_related_samples_to_drop.ht()
    relatedness_ht = get_relatedness_set_ht(hl.read_table(os.path.join(BASE, WRITE_DIR, "relatedness_annotated.ht"))) # relatedness.ht()
    related_samples_to_drop_ht = related_samples_to_drop_ht.annotate(
        relationships=relatedness_ht[related_samples_to_drop_ht.s].relationships
    )
    release_related_samples_to_drop_ht = release_related_samples_to_drop_ht.annotate(
        relationships=relatedness_ht[release_related_samples_to_drop_ht.s].relationships
    )

    # Annotating meta HT with related filter booleans
    # Any sample that is hard filtered will have missing values for these bools
    # Any sample that was filtered for relatedness will have True for sample_filters.related
    # If a filtered related sample had a relationship with a higher degree than second-degree (duplicate, parent-child, sibling),
    # that filter will also be True
    release_else_expr = release_related_samples_to_drop_ht[left_ht.key].relationships
    all_else_expr = related_samples_to_drop_ht[left_ht.key].relationships
    left_ht = left_ht.annotate(
        sample_filters=left_ht.sample_filters.annotate(
            release_related=hl.if_else(
                left_ht.sample_filters.hard_filtered,
                hl.null(hl.tbool),
                hl.is_defined(release_related_samples_to_drop_ht[left_ht.key]),
            ),
            release_duplicate=get_relationship_filter_expr(
                left_ht.sample_filters.hard_filtered,
                DUPLICATE_OR_TWINS,
                release_else_expr,
            ),
            release_parent_child=get_relationship_filter_expr(
                left_ht.sample_filters.hard_filtered,
                PARENT_CHILD,
                release_else_expr,
            ),
            release_sibling=get_relationship_filter_expr(
                left_ht.sample_filters.hard_filtered,
                SIBLINGS,
                release_else_expr,
            ),
            all_samples_related=hl.if_else(
                left_ht.sample_filters.hard_filtered,
                hl.null(hl.tbool),
                hl.is_defined(related_samples_to_drop_ht[left_ht.key]),
            ),
            all_samples_duplicate=get_relationship_filter_expr(
                left_ht.sample_filters.hard_filtered,
                DUPLICATE_OR_TWINS,
                all_else_expr,
            ),
            all_samples_parent_child=get_relationship_filter_expr(
                left_ht.sample_filters.hard_filtered,
                PARENT_CHILD,
                all_else_expr,
            ),
            all_samples_sibling=get_relationship_filter_expr(
                left_ht.sample_filters.hard_filtered,
                SIBLINGS,
                all_else_expr,
            ),
        )
    )
    left_ht = left_ht.annotate(
        relatedness_inference=hl.struct(
            relationships=relatedness_ht[left_ht.s].relationships,
        )
    )

    logger.info("Adding relatedness globals (cutoffs)")
    left_ht = left_ht.annotate_globals(
        relatedness_inference_cutoffs=hl.struct(**relatedness_ht.index_globals())
    )
    left_ht = left_ht.checkpoint(os.path.join(BASE, WRITE_DIR, "relatedness_joined.ht"), overwrite=True)

    logger.info(logging_statement.format("outlier HT"))
    if regressed_metrics_outlier:
        right_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, "regressed_metrics.ht")) # regressed_metrics.ht()
    else:
        right_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, "stratified_metrics.ht")) # stratified_metrics.ht()

    right_ht = right_ht.select_globals(
        outlier_detection_metrics=hl.struct(
            **right_ht.index_globals(), used_regressed_metrics=regressed_metrics_outlier
        )
    )

    left_ht = join_tables(left_ht, "s", right_ht, "s", "outer")
    left_ht = left_ht.transmute(
        sample_filters=left_ht.sample_filters.annotate(
            **{x: left_ht[x] for x in left_ht.row if x.startswith("fail_")},
            qc_metrics_filters=left_ht.qc_metrics_filters,
        )
    )
    if regressed_metrics_outlier:
        left_ht = left_ht.transmute(
            sample_qc=left_ht.sample_qc.annotate(
                **{x: left_ht[x] for x in left_ht.row if x.endswith("_residual")},
            )
        )
    left_ht = left_ht.checkpoint(os.path.join(BASE, WRITE_DIR, "outlier_joined.ht"), overwrite=True)

    logger.info("Annotating high_quality field")
    left_ht = left_ht.annotate(
        high_quality=~left_ht.sample_filters.hard_filtered
        & (hl.len(left_ht.sample_filters.qc_metrics_filters) == 0)
    )

    logger.info("Annotating releasable field")
    left_ht = left_ht.annotate(
        release=left_ht.project_meta.releasable
        & left_ht.high_quality
        & ~left_ht.project_meta.exclude
        & ~left_ht.sample_filters.release_related
    ).persist()

    logger.info(
        "Release sample counts:"
        f"{left_ht.aggregate(hl.struct(release=hl.agg.count_where(left_ht.release)))}"
    )
    left_ht.describe()
    #left_ht.summarize()
    logger.info(f"Final count: {left_ht.count()}")
    logger.info("Complete")

    return left_ht


def main(args):  # noqa: D103

    if args.sample_qc:
        compute_sample_qc().write(os.path.join(BASE, WRITE_DIR, 'sample_qc.ht'), overwrite=True)

    if args.impute_sex:
        compute_sex().write(os.path.join(BASE, WRITE_DIR, 'impute_sex.ht'), overwrite=True)
    elif args.reannotate_sex:
        reannotate_sex(
            args.min_cov,
            (args.upper_x, (args.lower_xx, args.upper_xx), args.lower_xxx),
            ((args.lower_y, args.upper_y), args.lower_yy),
        ).write(os.path.join(BASE, WRITE_DIR, 'reannotate_sex.ht'), overwrite=True)

    if args.compute_hard_filters:
        compute_hard_filters(
            args.min_cov,
            max_n_snp=args.max_n_snp,
            min_n_snp=args.min_n_snp,
            max_n_singleton=args.max_n_singleton,
            max_r_het_hom_var=args.max_r_het_hom_var,
        ).write(
            os.path.join(BASE, WRITE_DIR, 'compute_hard_filters.ht'), overwrite=True
        )

    if args.compute_qc_mt:
        compute_qc_mt().write(os.path.join(BASE, WRITE_DIR, 'compute_qc_mt.mt'), overwrite=True)

    if args.run_pc_relate or args.reannotate_relatedness:
        if args.run_pc_relate:
            logger.info("Running PC-Relate")
            logger.warning(
                "PC-relate requires SSDs and doesn't work with preemptible workers!"
            )
            # MT formerly being fed to hwe_normalized_pca (~50,000,000 variants)
            qc_mt = hl.read_matrix_table(os.path.join(BASE, WRITE_DIR, 'compute_qc_mt.mt'))
            # Keep biallelic sites only
            qc_mt = qc_mt.filter_rows(hl.len(qc_mt.alleles) == 2)
            # Compute variant_qc and filter by MAF >= 0.01
            qc_mt = hl.variant_qc(qc_mt)
            qc_mt = qc_mt.filter_rows(qc_mt.variant_qc.AF[0] >= 0.01)
            qc_mt = qc_mt.checkpoint(os.path.join(BASE, WRITE_DIR, "qc_mt_filtered.mt"), overwrite=True)
            # Take 0.1-1% of variants to downsample from ~50,000,000 to ~50,000-500,000
            sub_sample_mt = qc_mt.sample_rows(0.002, seed=42)
            # Selecting the relevant rows:
            selected_rows = sub_sample_mt.rows() #.select('locus', 'alleles')
            # Write sub-sampled table to disk
            selected_rows.write(os.path.join(BASE, WRITE_DIR, 'sub_sampled_qc.ht'), overwrite=True)
            # Read sub-sampled table from disk
            sub_sample_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, 'sub_sampled_qc.ht'))
            # Use semi_join_rows to get filtered mt to feed to hwe_normalized_pca
            filtered_qc_mt = qc_mt.semi_join_rows(sub_sample_ht)

            eig, scores, _ = hl.hwe_normalized_pca(
                filtered_qc_mt.GT, k=10, compute_loadings=False
            )
            scores = scores.checkpoint(
                os.path.join(BASE, WRITE_DIR, 'qc_mt_sites_pc_scores.ht'),
                overwrite=args.overwrite,
                _read_if_exists=not args.overwrite,
            )
            # Set min_kinship higher or lower to get more or less related samples (default 0.05) and/or min_individual_maf (default 0.01) 
            relatedness_ht = hl.pc_relate(
                filtered_qc_mt.GT,
                min_individual_maf=0.01,
                scores_expr=scores[filtered_qc_mt.col_key].scores,
                block_size=2048,
                min_kinship=0.05,
                statistics="all",
            )

        else:
            relatedness_ht = relatedness.ht().checkpoint(
                os.path.join(BASE, WRITE_DIR, 'run_pc_relate_checkpoint.ht'), overwrite=True
            )  # Copy HT to temp location to overwrite annotation
        relatedness_ht = relatedness_ht.annotate(
            relationship=get_relationship_expr(
                kin_expr=relatedness_ht.kin,
                ibd0_expr=relatedness_ht.ibd0,
                ibd1_expr=relatedness_ht.ibd1,
                ibd2_expr=relatedness_ht.ibd2,
                first_degree_kin_thresholds=tuple(args.first_degree_kin_thresholds),
                second_degree_min_kin=args.second_degree_kin_cutoff,
                ibd0_0_max=args.ibd0_0_max,
            )
        )
        relatedness_ht = relatedness_ht.annotate_globals(
            min_individual_maf=0.01,
            min_emission_kinship=0.05,
            ibd0_0_max=args.ibd0_0_max,
            second_degree_kin_cutoff=args.second_degree_kin_cutoff,
            first_degree_kin_thresholds=tuple(args.first_degree_kin_thresholds),
        )
        relatedness_ht.write(os.path.join(BASE, WRITE_DIR, 'run_pc_relate.ht'), True)

    if args.run_pca:
        rank_ht = compute_sample_rankings(
            use_qc_metrics_filters=False
        )  # QC metrics filters do not exist at this point
        rank_ht = rank_ht.checkpoint(
            os.path.join(BASE, WRITE_DIR, 'run_pca_rank.ht'),
            overwrite=args.overwrite,
            _read_if_exists=not args.overwrite,
        )
        filtered_samples = hl.literal(
            rank_ht.aggregate(
                hl.agg.filter(rank_ht.filtered, hl.agg.collect_as_set(rank_ht.s))
            ),
            dtype=hl.tset(hl.tstr)
        )  # TODO: don't localize once hail bug is fixed
        relatedness_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, 'run_pc_relate.ht'))
        relatedness_ht = relatedness_ht.key_by(
            i=relatedness_ht.key[0].s,
            j=relatedness_ht.key[1].s
        )
        samples_to_drop = compute_related_samples_to_drop(
            relatedness_ht,
            rank_ht,
            args.second_degree_kin_cutoff,
            filtered_samples=filtered_samples,
        )
        samples_to_drop = samples_to_drop.key_by(s=samples_to_drop.s)
        samples_to_drop.checkpoint(
            os.path.join(BASE, WRITE_DIR, 'samples_to_drop.ht'),
            overwrite=args.overwrite,
            _read_if_exists=not args.overwrite,
        )
        pop_pca_eigenvalues, pop_pca_scores_ht, pop_pca_loadings_ht = run_pca(
            args.include_unreleasable_samples, 10, samples_to_drop
        )
        pop_pca_scores_ht.write(
            os.path.join(BASE, WRITE_DIR, 'ancestry_pca_scores.ht'),
            overwrite=True,
        )
        pop_pca_loadings_ht.write(
            os.path.join(BASE, WRITE_DIR, 'ancestry_pca_loadings.ht'),
            overwrite=True,
        )

        pop_pca_eigenvalues_ht = hl.Table.parallelize(
            hl.literal(
                [
                    {"PC": i + 1, "eigenvalue": x}
                    for i, x in enumerate(pop_pca_eigenvalues)
                ],
                "array<struct{PC: int, eigenvalue: float}>",
            )
        )
        pop_pca_eigenvalues_ht.write(
            os.path.join(BASE, WRITE_DIR, 'pop_pca_eigenvalues.ht'),
            overwrite=True,
        )

    if args.assign_pops:
        pcs = args.pop_pcs
        pop_ht, pops_rf_model = assign_pops(
            args.min_pop_prob,
            args.include_unreleasable_samples,
            max_number_mislabeled_training_samples = 100,
            pcs=pcs,
            withhold_prop=args.withhold_prop,
        )
        pop_ht = pop_ht.checkpoint(
            os.path.join(BASE, WRITE_DIR, "pop.ht"),  # pop.path
            overwrite=args.overwrite,
            _read_if_exists=not args.overwrite,
        )
        # pca_scores holds only the PCs in `pcs`, in order, so index by position.
        pop_ht.transmute(
            **{f"PC{pc}": pop_ht.pca_scores[idx] for idx, pc in enumerate(pcs)}
        ).export(
            os.path.join(BASE, WRITE_DIR, "RF_pop_assignments.txt.gz")  # pop_tsv_path()
        )

        with hl.hadoop_open(
            os.path.join(BASE, WRITE_DIR, "pop.RF_fit.pickle"), "wb"  # pop_rf_path()
        ) as out:
            pickle.dump(pops_rf_model, out)

    if args.calculate_inbreeding:
        qc_mt = qc.mt()
        pop_ht = pop.ht()
        qc_mt = qc_mt.annotate_cols(pop=pop_ht[qc_mt.col_key].pop)
        qc_mt = qc_mt.annotate_rows(
            call_stats_by_pop=hl.agg.group_by(
                qc_mt.pop, hl.agg.call_stats(qc_mt.GT, qc_mt.alleles)
            )
        )
        inbreeding_ht = (
            qc_mt.annotate_cols(
                inbreeding=hl.agg.inbreeding(
                    qc_mt.GT, qc_mt.call_stats_by_pop[qc_mt.pop].AF[1]
                )
            )
            .cols()
            .select("inbreeding")
        )
        inbreeding_ht.write(sample_inbreeding.path, overwrite=args.overwrite)

    if args.calculate_clinvar:
        mt = get_gnomad_v3_mt(
            split=True,
            key_by_locus_and_alleles=True,
            remove_hard_filtered_samples=True
        )
        clinvar_ht = hl.read_table(os.path.join(BASE, 'hail', 'clinvar_20190923.ht'))

        # Normalize to minimal representation (left-aligned and trimmed)
        mt = mt.annotate_rows(alleles_minrep = hl.min_rep(mt.locus, mt.alleles))
        mt = mt.key_rows_by(mt.locus, mt.alleles_minrep)

        clinvar_ht = clinvar_ht.annotate(alleles_minrep = hl.min_rep(clinvar_ht.locus, clinvar_ht.alleles))
        clinvar_ht = clinvar_ht.key_by(clinvar_ht.locus, clinvar_ht.alleles_minrep)

        mt = mt.filter_rows(hl.is_defined(clinvar_ht[mt.locus, mt.alleles_minrep]))
        mt = mt.checkpoint(
            os.path.join(BASE, WRITE_DIR, 'clinvar_variants.mt'),
            overwrite=True
        )

        path_ht = filter_to_clinvar_pathogenic(clinvar_ht)
        path_ht = path_ht.key_by(
            locus=path_ht.locus,
            alleles_minrep=hl.min_rep(path_ht.locus, path_ht.alleles)
        )

        mt = mt.annotate_rows(
            clinvar_path = hl.is_defined(path_ht[mt.locus, mt.alleles_minrep])
        )

        # --- Aggregate counts per sample ---
        clinvar_sample_ht = mt.select_cols(
            n_clinvar=hl.agg.count_where(mt.GT.is_non_ref()),
            n_clinvar_path=hl.agg.count_where(mt.GT.is_non_ref() & mt.clinvar_path),
        ).cols()

        # Write per-sample ClinVar counts
        clinvar_sample_ht.write(
            os.path.join(BASE, WRITE_DIR, 'sample_clinvar_count.ht'),
            overwrite=args.overwrite
        )

        # --- (Optional) Sanity checks ---
        n_rows_total = mt.count_rows()
        n_path_rows = mt.aggregate_rows(hl.agg.count_where(mt.clinvar_path))
        print(f"Total ClinVar-overlapping variants: {n_rows_total}")
        print(f"Pathogenic ClinVar variants: {n_path_rows}")

    if args.apply_stratified_filters or args.apply_regressed_filters:
        filtering_qc_metrics = args.filtering_qc_metrics.split(",")
        sample_qc_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, "tmp_sample_qc.ht_bi_allelic.ht")) # get_sample_qc("bi_allelic").ht()

        if "inbreeding" in filtering_qc_metrics:
            inbreeding_ht = sample_inbreeding.ht()[sample_qc_ht.key]
            sample_qc_ht = sample_qc_ht.annotate(
                sample_qc=sample_qc_ht.sample_qc.annotate(
                    inbreeding=inbreeding_ht.inbreeding.f_stat
                )
            )

        if args.apply_regressed_filters:
            n_pcs = args.regress_n_pcs
            apply_regressed_filters(
                sample_qc_ht,
                filtering_qc_metrics,
                args.include_unreleasable_samples,
                n_pcs,
            ).write(os.path.join(BASE, WRITE_DIR, "regressed_metrics.ht"), overwrite=args.overwrite)
        else:
            apply_stratified_filters(
                sample_qc_ht,
                filtering_qc_metrics,
            ).write(os.path.join(BASE, WRITE_DIR, "stratified_metrics.ht"), overwrite=args.overwrite)

    if args.compute_related_samples_to_drop:
        rank_ht = compute_sample_rankings(use_qc_metrics_filters=True)
        rank_ht = rank_ht.checkpoint(
            os.path.join(BASE, WRITE_DIR, 'release_samples_rankings.ht'),
            overwrite=args.overwrite,
            _read_if_exists=not args.overwrite,
        )
        filtered_samples = hl.literal(
            rank_ht.aggregate(
                hl.agg.filter(rank_ht.filtered, hl.agg.collect_as_set(rank_ht.s))
            ),
            dtype=hl.tset(hl.tstr)
        )  # TODO: don't localize once hail bug is fixed
        relatedness_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, 'run_pc_relate.ht'))
        relatedness_ht = relatedness_ht.key_by(
            i=relatedness_ht.i.s, j=relatedness_ht.j.s
        )
        samples_to_drop = compute_related_samples_to_drop(
            relatedness_ht,
            rank_ht,
            args.second_degree_kin_cutoff,
            filtered_samples=None,
        )
        samples_to_drop.write(
            os.path.join(BASE, WRITE_DIR, 'release_related_samples_to_drop.ht'), overwrite=args.overwrite
        )

    if args.generate_metadata:
        meta_ht = generate_metadata(args.regressed_metrics_outlier)
        final_meta = os.path.join(BASE, WRITE_DIR, 'final_meta')
        meta_ht.checkpoint(
            final_meta + ".ht", overwrite=args.overwrite, _read_if_exists=not args.overwrite
        )
        n_pcs = meta_ht.aggregate(
            hl.agg.min(hl.len(meta_ht.population_inference.pca_scores))
        )

        meta_ht = meta_ht.annotate(
            population_inference=meta_ht.population_inference.annotate(
                **{
                    f"PC{i + 1}": meta_ht.population_inference.pca_scores[i]
                    for i in range(n_pcs)
                }
            ),
            hard_filters=hl.or_missing(
                hl.len(meta_ht.sample_filters.hard_filters) > 0,
                hl.delimit(meta_ht.sample_filters.hard_filters),
            ),
            qc_metrics_filters=hl.or_missing(
                hl.len(meta_ht.sample_filters.qc_metrics_filters) > 0,
                hl.delimit(meta_ht.sample_filters.qc_metrics_filters),
            ),
        )
        meta_ht.flatten().export(final_meta + ".tsv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--sample_qc", help="Assigns pops from PCA", action="store_true"
    )
    parser.add_argument(
        "--impute_sex",
        help="Runs sex imputation. Also runs sex karyotyping annotation.",
        action="store_true",
    )
    parser.add_argument(
        "--reannotate_sex",
        help=(
            "Runs the sex karyotyping annotations again, without re-computing sex"
            " imputation metrics."
        ),
        action="store_true",
    )
    parser.add_argument("--upper_x", help="Upper cutoff for single X", type=float)
    parser.add_argument("--lower_xx", help="Lower cutoff for double X", type=float)
    parser.add_argument("--upper_xx", help="Upper cutoff for double X", type=float)
    parser.add_argument("--lower_xxx", help="Lower cutoff for triple X", type=float)
    parser.add_argument("--lower_y", help="Lower cutoff for single Y", type=float)
    parser.add_argument("--upper_y", help="Upper cutoff for single Y", type=float)
    parser.add_argument("--lower_yy", help="Lower cutoff for double Y", type=float)
    parser.add_argument(
        "--compute_hard_filters",
        help="Computes samples to be hard-filtered",
        action="store_true",
    )
    parser.add_argument(
        "--min_cov",
        help="Minimum coverage for inclusion when computing hard-filters",
        default=10,
        type=int,
    )
    # Hard-filter cutoffs. The defaults are genome-wide counts; lower --min_n_snp
    # for an exome, panel or subset, or every sample is filtered out.
    parser.add_argument(
        "--max_n_snp",
        help="Maximum number of SNPs for a sample to pass hard filters",
        default=2e6,
        type=float,
    )
    parser.add_argument(
        "--min_n_snp",
        help="Minimum number of SNPs for a sample to pass hard filters",
        default=5e5,
        type=float,
    )
    parser.add_argument(
        "--max_n_singleton",
        help="Maximum number of singletons for a sample to pass hard filters",
        default=5e5,
        type=float,
    )
    parser.add_argument(
        "--max_r_het_hom_var",
        help="Maximum het/hom-var ratio for a sample to pass hard filters",
        default=11,
        type=float,
    )
    parser.add_argument(
        "--compute_qc_mt",
        help="Creates the QC MT based on liftover of v2 QC and Purcell 5k sites",
        action="store_true",
    )
    parser.add_argument("--run_pc_relate", help="Run PC-relate", action="store_true")
    parser.add_argument(
        "--reannotate_relatedness",
        help="Runs the relatedness annotation without re-running pc-relate",
        action="store_true",
    )
    parser.add_argument(
        "--first_degree_kin_thresholds",
        help=(
            "First degree kinship threshold for filtering a pair of samples with a"
            " first degree relationship.         Default = (0.1767767, 0.4);        "
            " Defaults taken from Bycroft et al. (2018)"
        ),
        nargs=2,
        default=(0.1767767, 0.4),
        type=float,
    )
    parser.add_argument(
        "--second_degree_kin_cutoff",
        help=(
            "Minimum kinship threshold for filtering a pair of samples with a second"
            " degree relationship        in PC relate and filtering related"
            " individuals. (Default = 0.1)         Bycroft et al. (2018) calculates"
            " 0.08838835 but from evaluation of the distributions v3 has used 0.1"
        ),
        default=0.1,
        type=float,
    )
    parser.add_argument(
        "--ibd0_0_max",
        help=(
            "IBD0 cutoff to determine parent offspring vs full sibling (Default = 0.05)"
            "         Default is adjusted from theoretical values; parent-offspring"
            " should have an IBD0 = 0.         Full siblings should have an IBD0 ="
            " 0.25."
        ),
        default=0.05,
    )
    parser.add_argument(
        "--compute_related_samples_to_drop",
        help="Flags related samples to drop",
        action="store_true",
    )
    parser.add_argument(
        "--min_related_hard_filter",
        help="Minimum number of relateds to have to get hard-filterd",
        default=50,
        type=int,
    )
    parser.add_argument("--run_pca", help="Compute PCA", action="store_true")
    parser.add_argument(
        "--n_pcs",
        help="Number of PCs to compute for ancestry PCA",
        default=30,
        type=int,
    )
    parser.add_argument(
        "--include_unreleasable_samples",
        help="Includes unreleasable samples for computing PCA",
        action="store_true",
    )
    parser.add_argument(
        "--assign_pops", help="Assigns pops from PCA", action="store_true"
    )
    parser.add_argument(
        "--pop_pcs",
        help=(
            "List of PCs to use for ancestry assignment. The values provided should be"
            " 1-based."
        ),
        default=list(range(1, 3)),
        # Space-separated, e.g. --pop_pcs 1 2 3
        type=int,
        nargs="+",
    )
    parser.add_argument(
        "--min_pop_prob",
        help="Minimum RF prob for pop assignment",
        default=0.75,
        type=float,
    )
    parser.add_argument(
        "--withhold_prop",
        help=(
            "Proportion of training samples to withhold from pop assignment RF training"
        ),
        type=float,
    )
    parser.add_argument(
        "--calculate_inbreeding",
        help="Calculate sample level inbreeding",
        action="store_true",
    )
    parser.add_argument(
        "--calculate_clinvar",
        help="Calculate counts of ClinVar and ClinVar P/LP variants per sample",
        action="store_true",
    )
    parser.add_argument(
        "--filtering_qc_metrics",
        help="List of QC metrics for filtering.",
        default=",".join(
            [
                "n_snp",
                "n_singleton",
                "r_ti_tv",
                "r_insertion_deletion",
                "n_insertion",
                "n_deletion",
                "r_het_hom_var",
                "n_transition",
                "n_transversion",
            ]
        ),
    )  # used in v3 'n_het', 'n_hom_var',
    parser.add_argument(
        "--apply_stratified_filters",
        help="Compute per pop filtering.",
        action="store_true",
    )
    parser.add_argument(
        "--apply_regressed_filters",
        help="Computes qc_metrics adjusted for pop.",
        action="store_true",
    )
    parser.add_argument(
        "--regress_n_pcs",
        help="Number of PCs to use for qc metric regressions",
        default=10,
        type=int,
    )
    parser.add_argument(
        "--generate_metadata", help="Generates the metadata HT.", action="store_true"
    )
    parser.add_argument(
        "--regressed_metrics_outlier",
        help="Should metadata HT use regression outlier model.",
        action="store_true",
    )
    parser.add_argument(
        "--seed",
        help="Random seed",
        default=24,
        type=int,
    )

    main(parser.parse_args())
# %%