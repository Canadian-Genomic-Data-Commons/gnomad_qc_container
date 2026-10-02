# noqa: D100

import argparse
import logging
from collections import Counter, defaultdict

import os
import hail as hl

from gnomad_qc.pipeline_config import BASE, WRITE_DIR, init_hail
from gnomad_qc.v3.resources.basics import get_gnomad_v3_mt
from gnomad_qc.v3.resources.constants import CURRENT_VERSION
from gnomad_qc.v3.resources.meta import trios
from gnomad_qc.v3.resources.sample_qc import (
    duplicates,
    get_relatedness_annotated_ht,
    release_samples_rankings,
    sex,
)
from gnomad_qc.v3.resources.relatedness import (
    create_fake_pedigree,
    get_duplicated_samples,
    get_duplicated_samples_ht,
    infer_families,
)
from gnomad.sample_qc.relatedness import get_relationship_expr
logger = logging.getLogger("create_fam")


def run_mendel_errors() -> hl.Table:  # noqa: D103
    meta_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, "meta.ht"))
    ped = hl.Pedigree.read(os.path.join(BASE, WRITE_DIR, "pedigree_raw.fam")) #pedigree.versions[f"{CURRENT_RELEASE}_raw"].pedigree()
    logger.info(f"Running Mendel errors for {len(ped.trios)} trios.")
    """
    fake_ped = create_fake_pedigree(
        n=100,
        sample_list=list(
            meta_ht.aggregate(
                hl.agg.filter(
                    hl.rand_bool(0.01)
                    & (
                        (hl.len(meta_ht.qc_metrics_filters) == 0)
                        & hl.or_else(hl.len(meta_ht.hard_filters) == 0, False)
                    ),
                    hl.agg.collect_as_set(meta_ht.s),
                )
            )
        ),
        real_pedigree=ped,
    )
    """
    fake_ped = create_fake_pedigree(
        n=100,
        sample_list=list(
            meta_ht.aggregate(
                hl.agg.filter(
                    hl.rand_bool(0.01),
                    hl.agg.collect_as_set(meta_ht.s),
                )
            )
        ),
        real_pedigree=ped,
    )
    merged_ped = hl.Pedigree(trios=ped.trios + fake_ped.trios)

    ped_samples = hl.literal(
        set(
            [s for trio in merged_ped.trios for s in [trio.s, trio.pat_id, trio.mat_id]]
        )
    )
    mt = get_gnomad_v3_mt(key_by_locus_and_alleles=True)
    mt = mt.filter_cols(ped_samples.contains(mt.s))
    mt = hl.filter_intervals(
        mt, [hl.parse_locus_interval("chr20", reference_genome="GRCh38")]
    )
    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)
    mt = mt.select_entries("GT", "END")
    mt = hl.experimental.densify(mt)
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)
    mendel_errors, _, _, _ = hl.mendel_errors(mt["GT"], merged_ped)
    return mendel_errors


def run_infer_families() -> hl.Pedigree:
    """
    Note: if there are few rows in run_pc_relate.ht table there may be no trios inferred. 
    Change settings such as min_kinship and the degree of subsampling in the pc_relate function in the relatedness.py script to get more trios.
    """
    if os.path.exists(os.path.join(BASE, WRITE_DIR, 'duplicates.ht')):
        dups_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, 'duplicates.ht'))
    else:
        def mock_empty_duplicated_samples_ht():
            return hl.Table.parallelize(
                [],
                schema=hl.tstruct(s=hl.tstr, filtered=hl.tarray(hl.tstr)),
                key='s'
            )
        dups_ht = mock_empty_duplicated_samples_ht()

    ped = infer_families(
        hl.read_table(os.path.join(BASE, WRITE_DIR, 'relatedness_annotated.ht')),
        hl.read_table(os.path.join(BASE, WRITE_DIR, 'impute_sex.ht')),
        dups_ht
    )
    print(f"Initial trios inferred: {len(ped.trios)}")

    meta_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, "meta.ht"))
    filtered_samples = meta_ht.aggregate(
        hl.agg.filter(
            hl.is_missing(meta_ht.group_by_expr),
            hl.agg.collect_as_set(meta_ht.s),
        )
    )
    print(f"Filtered samples count: {len(filtered_samples)}")
    print(f"Filtered samples: {filtered_samples}")

    filtered_trios = [
        trio
        for trio in ped.trios
        if trio.s not in filtered_samples
        and trio.pat_id not in filtered_samples
        and trio.mat_id not in filtered_samples
    ]
    print(f"Number of trios retained after filtering: {len(filtered_trios)}")
    return hl.Pedigree(trios=filtered_trios)



def families_to_trios(ped: hl.Pedigree) -> hl.Pedigree:  # noqa: D103
    trios_per_fam = defaultdict(list)
    for trio in ped.trios:
        trios_per_fam[trio.fam_id].append(trio)

    message = [
        f"Found a total of {len(ped.trios)} trios in {len(trios_per_fam)} families:"
    ]
    for n_trios, n_fam in Counter(
        [len(trios) for trios in trios_per_fam.values()]
    ).items():
        message.append(f"{n_fam} with {n_trios} trios.")
    logger.info("\n".join(message))

    return hl.Pedigree(trios=[trios[0] for trios in trios_per_fam.values()])


def filter_ped(
    raw_ped: hl.Pedigree, mendel: hl.Table, max_dnm: int, max_mendel: int
) -> hl.Pedigree:  # noqa: D103
    mendel = mendel.filter(mendel.fam_id.startswith("fake"))
    mendel_by_s = (
        mendel.group_by(mendel.s)
        .aggregate(
            fam_id=hl.agg.take(mendel.fam_id, 1)[0],
            n_mendel=hl.agg.count(),
            n_de_novo=hl.agg.count_where(
                mendel.mendel_code == 2
            ),  # Code 2 is parents are hom ref, child is het
        )
        .persist()
    )

    good_trios = mendel_by_s.aggregate(
        hl.agg.filter(
            (mendel_by_s.n_mendel < max_mendel) & (mendel_by_s.n_de_novo < max_dnm),
            hl.agg.collect(
                mendel_by_s.s,
            ),
        )
    )
    logger.info(f"Found {len(good_trios)} trios passing filters")
    return hl.Pedigree([trio for trio in raw_ped.trios if trio.s in good_trios])


def main(args):  # noqa: D103
    init_hail("create_fam")

    if args.find_dups:
        logger.info("Selecting best duplicate per duplicated sample set")
        relatedness_ht = hl.read_table(os.path.join(BASE, WRITE_DIR, 'run_pc_relate.ht')) # get_relatedness_annotated_ht()
        relatedness_annotated = relatedness_ht.annotate(
        relationship=get_relationship_expr(
            kin_expr=relatedness_ht.kin,
            ibd0_expr=relatedness_ht.ibd0,
            ibd1_expr=relatedness_ht.ibd1,
            ibd2_expr=relatedness_ht.ibd2,
            )
        )
        relatedness_annotated.write(os.path.join(BASE, WRITE_DIR, 'relatedness_annotated.ht'), overwrite=args.overwrite)
        dups = get_duplicated_samples(relatedness_annotated)
        dups_ht = get_duplicated_samples_ht(dups, hl.read_table(os.path.join(BASE, WRITE_DIR, 'release_samples_rankings.ht'))) # release_samples_rankings.ht()
        if dups_ht.count() == 0:
            logger.warning("The duplicates HT is empty. Skipping write operation.")
        else:
            dups_ht.write(os.path.join(BASE, WRITE_DIR, 'duplicates.ht'), overwrite=args.overwrite)
    if args.infer_families:
        ped = run_infer_families()
        ped.write(os.path.join(BASE, WRITE_DIR, 'pedigree_raw.fam'))
        raw_trios = families_to_trios(ped)
        raw_trios.write(os.path.join(BASE, WRITE_DIR, 'trios_raw.fam'))

    if args.run_mendel_errors:
        mendel_errors = run_mendel_errors()
        # Table.write takes a path. This was passing hl.Pedigree.read() of the
        # output path, so the step tried to read the file it was about to write
        # and always failed with FileNotFoundException on a first run.
        mendel_errors.write(
            os.path.join(BASE, WRITE_DIR, "ped_mendel_errors.ht"),
            overwrite=args.overwrite,
        )

    if args.finalize_ped:
        final_ped = filter_ped(
            hl.Pedigree.read(os.path.join(BASE, WRITE_DIR, "pedigree_raw.fam")),
            hl.read_table(os.path.join(BASE, WRITE_DIR, "ped_mendel_errors.ht")), # ped_mendel_errors.ht()
            args.max_dnm,
            args.max_mendel,
        )
        final_ped.write(os.path.join(BASE, WRITE_DIR, 'pedigree_final.fam'))

        final_trios = families_to_trios(final_ped)
        final_trios.write(trios.path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--find_dups",
        help=(
            "Creates a table with duplicate samples indicating which one is the best to"
            " use."
        ),
        action="store_true",
    )
    parser.add_argument("--infer_families", help="Infers families", action="store_true")
    parser.add_argument(
        "--run_mendel_errors", help="Runs mendel errors", action="store_true"
    )
    parser.add_argument(
        "--finalize_ped",
        help=(
            "Creates a final ped file by excluding families where the number of Mendel"
            " errors or de novos are higher than those specified in --max_dnm and"
            " --max_mendel"
        ),
        action="store_true",
    )
    parser.add_argument(
        "--max_dnm",
        help="Maximum number of raw de novo mutations for real trios",
        default=2200,
        type=int,
    )
    parser.add_argument(
        "--max_mendel",
        help="Maximum number of raw Mendel errors for real trios",
        default=3750,
        type=int,
    )

    main(parser.parse_args())
