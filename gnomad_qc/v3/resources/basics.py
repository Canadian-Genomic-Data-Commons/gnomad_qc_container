# noqa: D100
import logging
from typing import List, Optional, Set, Union
import os
import hail as hl
from gnomad.resources.resource_utils import (
    MatrixTableResource,
    VariantDatasetResource,
    VersionedMatrixTableResource,
    VersionedVariantDatasetResource,
)

from gnomad_qc.pipeline_config import BASE, WRITE_DIR, VDS_FILE

#import gnomad_qc.v4.resources.basics as v4_basics
from gnomad_qc.v3.resources.constants import CURRENT_RELEASE, CURRENT_VERSION
from gnomad_qc.v3.resources.meta import meta
from gnomad_qc.v3.resources.sample_qc import hard_filtered_samples

logger = logging.getLogger("basic_resources")
logger.setLevel(logging.INFO)

def get_gnomad_v3_vds(
    split: bool = False,
    remove_hard_filtered_samples: bool = True,
    release_only: bool = False,
    samples_meta: bool = False,
    test: bool = False,
    filter_partitions: Optional[List[int]] = None,
    chrom: Optional[Union[str, List[str], Set[str]]] = None,
    autosomes_only: bool = False,
    sex_chr_only: bool = False,
    filter_variant_ht: Optional[hl.Table] = None,
    filter_intervals: Optional[List[Union[str, hl.tinterval]]] = None,
    split_reference_blocks: bool = True,
    entries_to_keep: Optional[List[str]] = None,
    annotate_het_non_ref: bool = False,
    naive_coalesce_partitions: Optional[int] = None,
    filter_samples_ht: Optional[hl.Table] = None,
) -> hl.vds.VariantDataset:
    """
    Get gnomAD VariantDataset with desired filtering and metadata annotations.

    :param split: Perform split on VDS - Note: this will perform a split on the VDS
        rather than grab an already split VDS.
    :param remove_hard_filtered_samples: Whether to remove samples that failed hard
        filters (only relevant after sample QC).
    :param release_only: Whether to filter the VDS to only samples available for
        release (can only be used if metadata is present).
    :param samples_meta: Whether to add metadata to VDS variant_data in 'meta' column.
    :param test: Whether to use the test VDS instead of the full v3 VDS.
    :param filter_partitions: Optional argument to filter the VDS to specific partitions
        in the provided list.
    :param chrom: Optional argument to filter the VDS to specific chromosomes.
    :param autosomes_only: Whether to filter the VDS to autosomes only. Default is
        False.
    :param sex_chr_only: Whether to filter the VDS to sex chromosomes only. Default is
        False.
    :param filter_variant_ht: Optional argument to filter the VDS to a specific set of
        variants. Only supported when splitting the VDS.
    :param filter_intervals: Optional argument to filter the VDS to specific intervals.
    :param split_reference_blocks: Whether to split the reference data at the edges of
        the intervals defined by `filter_intervals`. Default is True.
    :param entries_to_keep: Optional argument to keep only specific entries in the
        returned VDS. If splitting the VDS, use the global entries (e.g. 'GT') instead
        of the local entries (e.g. 'LGT') to keep.
    :param annotate_het_non_ref: Whether to annotate non reference heterozygotes (as
        '_het_non_ref') to the variant data. Default is False.
    :param naive_coalesce_partitions: Optional argument to coalesce the VDS to a
        specific number of partitions using naive coalesce.
    :param filter_samples_ht: Optional Table of samples to filter the VDS to.
    :return: gnomAD v3 dataset with chosen annotations and filters.
    """
    vds=hl.vds.read_vds(os.path.join(BASE, VDS_FILE))

    if remove_hard_filtered_samples:
        vds = hl.vds.filter_samples(
            vds,
            hl.read_table(os.path.join(BASE, WRITE_DIR, 'compute_hard_filters.ht')),
            keep=False,
        )

    if filter_samples_ht:
        logger.info(
            "Filtering VDS to %d samples in provided Table...",
            filter_samples_ht.count(),
        )
        vds = hl.vds.filter_samples(vds, filter_samples_ht)

    if samples_meta or release_only:
        meta_ht = hl.read_table(meta)
        if release_only:
            vds = hl.vds.filter_samples(
                vds,
                meta_ht.filter(meta_ht.release),
            )

        if samples_meta:
            vd = vds.variant_data
            vd.describe()
            print(f"This is vd.col_key: {vd.col_key}")
            meta_ht = meta_ht.key_by('s')
            meta_ht[vd.col_key]
            vds = hl.vds.VariantDataset(
                vds.reference_data, vd.annotate_cols(meta=meta_ht[vd.col_key])
            )

    vd = vds.variant_data
    if annotate_het_non_ref:
        logger.info("Annotating non_ref hets to unsplit variant data...")
        vd = vd.annotate_entries(_het_non_ref=vd.LGT.is_het_non_ref())
        entries_to_keep = (
            None if entries_to_keep is None else entries_to_keep + ["_het_non_ref"]
        )

    #if split:
    #    vd = v4_basics._split_and_filter_variant_data_for_loading(
    #        vd, filter_variant_ht, entries_to_keep
    #    )

    if entries_to_keep is not None:
        vd = vd.select_entries(*entries_to_keep)

    vds = hl.vds.VariantDataset(vds.reference_data, vd)

    return vds


def get_gnomad_v3_mt(
    split: bool = False,
    key_by_locus_and_alleles: bool = False,
    remove_hard_filtered_samples: bool = True,
    release_only: bool = False,
    samples_meta: bool = False,
    test: bool = False,
) -> hl.MatrixTable:
    """
    Get gnomAD data with desired filtering and metadata annotations.

    :param split: Perform split on MT - Note: this will perform a split on the MT rather than grab an already split MT
    :param key_by_locus_and_alleles: Whether to key the MatrixTable by locus and alleles (only needed for v3)
    :param remove_hard_filtered_samples: Whether to remove samples that failed hard filters (only relevant after sample QC)
    :param release_only: Whether to filter the MT to only samples available for release (can only be used if metadata is present)
    :param samples_meta: Whether to add metadata to MT in 'meta' column
    :param test: Whether to use the test MT instead of the full v3 MT
    :return: gnomAD v3 dataset with chosen annotations and filters
    """
    logger.warning(
        "The raw v3.1 sparse MatrixTable has been deleted. This function now calls "
        "'get_gnomad_v3_vds' and turns the VariantDataset into a sparse MatrixTable."
        " We recommend direct use of the VariantDataset by calling 'get_gnomad_v3_vds' "
        "instead."
    )

    if test:
        mt = gnomad_v3_testset.mt()
        if key_by_locus_and_alleles:
            # Prevents hail from running sort on genotype MT which is already sorted
            # by a unique locus
            mt = hl.MatrixTable(
                hl.ir.MatrixKeyRowsBy(mt._mir, ["locus", "alleles"], is_sorted=True)
            )
        if remove_hard_filtered_samples:
            mt = mt.filter_cols(hl.is_missing(hard_filtered_samples.ht()[mt.col_key]))
        if samples_meta or release_only:
            meta_ht = hl.read_table(meta)
            if samples_meta:
                mt = mt.annotate_cols(meta=meta_ht[mt.col_key])
            if release_only:
                mt = mt.filter_cols(meta_ht[mt.col_key].release)
        if split:
            mt = mt.annotate_rows(
                n_unsplit_alleles=hl.len(mt.alleles),
                mixed_site=(hl.len(mt.alleles) > 2)
                & hl.any(lambda a: hl.is_indel(mt.alleles[0], a), mt.alleles[1:])
                & hl.any(lambda a: hl.is_snp(mt.alleles[0], a), mt.alleles[1:]),
            )
            mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)
    else:
        vds = get_gnomad_v3_vds(
            split=split,
            remove_hard_filtered_samples=remove_hard_filtered_samples,
            release_only=release_only,
            samples_meta=samples_meta,
        )
        mt = hl.vds.to_merged_sparse_mt(vds, ref_allele_function=lambda locus: hl.missing('str'))

        if not key_by_locus_and_alleles:
            mt = mt.key_rows_by(mt.locus)

    return mt


def qc_temp_prefix(version: str = CURRENT_RELEASE) -> str:
    """
    Return path to temporary QC bucket.

    :param version: Version of annotation path to return
    :return: Path to bucket with temporary QC data
    """
    return f"gs://gnomad-tmp-4day/gnomad_v{version}_qc_data/"


def get_checkpoint_path(
    name: str, version: str = CURRENT_RELEASE, mt: bool = False
) -> str:
    """
    Create a checkpoint path for Table or MatrixTable.

    :param str name: Name of intermediate Table/MatrixTable
    :param version: Version of annotation path to return
    :param bool mt: Whether path is for a MatrixTable, default is False
    :return: Output checkpoint path
    """
    return f'{qc_temp_prefix(version)}{name}.{"mt" if mt else "ht"}'


def get_logging_path(name: str, version: str = CURRENT_RELEASE) -> str:
    """
    Create a path for Hail log files.

    :param str name: Name of log file
    :param version: Version of annotation path to return
    :return: Output log path
    """
    return f"{qc_temp_prefix(version)}{name}.log"


# Resources need for SV histogram generation
gnomad_sv_bucket_path = "gs://talkowski-sv-gnomad-v3-release/Releasable_freeze1_202303"
"""
Path to bucket with gnomAD SV results (VCFs per chromosome).

Bucket also contains list of releasable samples (see below).
"""

gnomad_sv_release_samples_list_path = (
    f"{gnomad_sv_bucket_path}/gnomAD.v3.SV.releasable.samples"
)
