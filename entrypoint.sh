#!/usr/bin/env bash
# Entrypoint for the gnomad_qc container. Maps a --step name to the matching
# pipeline script (and the flags that step needs) and runs it.
#
#   apptainer run --bind /path/to/data:/data gnomad_qc.sif \
#       --step ancestry_classification
#
# Step names and their flags follow the pipeline order in the gnomAD SNV v3
# System Requirements and Benchmarking SOP. Extra arguments are passed through
# to the script, so e.g. `--step sample_qc --overwrite` works.
#
# GNOMAD_QC_CONFIG must point to a filled-in config.yaml (paths + Spark
# settings) that is bind-mounted into the container -- see
# gnomad_qc/v3/configuration/config.yaml for the template/schema.
set -euo pipefail

# Ordered for display; associative arrays below do not preserve order.
STEP_ORDER=(
    create_last_END_positions
    generate_qc_annotations
    generate_ac
    generate_fam_stats
    impute_sex
    sample_qc
    compute_hard_filters
    compute_qc_mt
    run_pca
    pc_relate
    assign_pops
    ancestry_classification
    generate_metadata
    find_dups
    infer_families
    run_mendel_errors
    finalize_ped
    annotate_for_rf
    train_rf
    apply_rf
    create_bin_ht
    create_aggregated_bin_ht
    final_filter
    freq_write_split_vds
    freq_run_dense_annotations
    freq_correct_high_ab_hets
    freq_finalize
    vep
    create_release_sites_ht
    generate_freq_data_v3
    make_pca_plot
    subpop_analysis
    v2_pc_relate
)

declare -A STEP_SCRIPT=(
    [create_last_END_positions]="gnomad_qc/v3/load_data/create_last_END_positions.py"
    [generate_qc_annotations]="gnomad_qc/v3/annotations/generate_qc_annotations.py"
    [generate_ac]="gnomad_qc/v3/annotations/generate_qc_annotations.py"
    [generate_fam_stats]="gnomad_qc/v3/annotations/generate_qc_annotations.py"
    [vep]="gnomad_qc/v3/annotations/generate_qc_annotations.py"
    [impute_sex]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [sample_qc]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [compute_hard_filters]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [compute_qc_mt]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [run_pca]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [pc_relate]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [assign_pops]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [generate_metadata]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [ancestry_classification]="gnomad_qc/v3/sample_qc/ancestry_classification.py"
    [find_dups]="gnomad_qc/v3/sample_qc/create_fam.py"
    [infer_families]="gnomad_qc/v3/sample_qc/create_fam.py"
    [run_mendel_errors]="gnomad_qc/v3/sample_qc/create_fam.py"
    [finalize_ped]="gnomad_qc/v3/sample_qc/create_fam.py"
    [annotate_for_rf]="gnomad_qc/v3/variant_qc/random_forest.py"
    [train_rf]="gnomad_qc/v3/variant_qc/random_forest.py"
    [apply_rf]="gnomad_qc/v3/variant_qc/random_forest.py"
    [create_bin_ht]="gnomad_qc/v3/variant_qc/evaluation.py"
    [create_aggregated_bin_ht]="gnomad_qc/v3/variant_qc/evaluation.py"
    [final_filter]="gnomad_qc/v3/variant_qc/final_filter.py"
    [freq_write_split_vds]="gnomad_qc/v4/annotations/generate_freq.py"
    [freq_run_dense_annotations]="gnomad_qc/v4/annotations/generate_freq.py"
    [freq_correct_high_ab_hets]="gnomad_qc/v4/annotations/generate_freq.py"
    [freq_finalize]="gnomad_qc/v4/annotations/generate_freq.py"
    [create_release_sites_ht]="gnomad_qc/v3/create_release/create_release_sites_ht.py"
    [generate_freq_data_v3]="gnomad_qc/v3/annotations/generate_freq_data.py"
    [make_pca_plot]="gnomad_qc/v3/sample_qc/make_pca_plot.py"
    [subpop_analysis]="gnomad_qc/v3/sample_qc/subpop_analysis.py"
    [v2_pc_relate]="gnomad_qc/v3/sample_qc/v2_pc_relate.py"
)

# Empty means the script takes no step-selecting flags (it runs top to bottom).
declare -A STEP_FLAGS=(
    [create_last_END_positions]=""
    [generate_qc_annotations]="--compute_info --split_info"
    [generate_ac]="--generate_ac"
    [generate_fam_stats]="--generate_fam_stats"
    [vep]="--vep"
    [impute_sex]="--impute_sex"
    [sample_qc]="--sample_qc"
    [compute_hard_filters]="--compute_hard_filters"
    [compute_qc_mt]="--compute_qc_mt"
    [run_pca]="--run_pca"
    [pc_relate]="--run_pc_relate"
    [assign_pops]="--assign_pops"
    [generate_metadata]="--generate_metadata"
    [ancestry_classification]=""
    [find_dups]="--find_dups"
    [infer_families]="--infer_families"
    [run_mendel_errors]="--run_mendel_errors"
    [finalize_ped]="--finalize_ped"
    [annotate_for_rf]="--annotate_for_rf"
    [train_rf]="--train_rf"
    [apply_rf]="--apply_rf"
    [create_bin_ht]="--create_bin_ht"
    [create_aggregated_bin_ht]="--create_aggregated_bin_ht"
    [final_filter]=""
    [freq_write_split_vds]="--write-split-vds-and-downsampling-ht"
    [freq_run_dense_annotations]="--run-freq-and-dense-annotations"
    [freq_correct_high_ab_hets]="--correct-for-high-ab-hets"
    [freq_finalize]="--finalize-freq-ht"
    [create_release_sites_ht]=""
    [generate_freq_data_v3]=""
    [make_pca_plot]=""
    [subpop_analysis]=""
    [v2_pc_relate]="--run_pc_relate"
)

usage() {
    cat <<EOF
Usage: --step <name> [extra args passed through to the script]

Steps, in SOP pipeline order (VDS creation happens before these, via build_vds/):

EOF
    for k in "${STEP_ORDER[@]}"; do
        printf '  %-28s %s %s\n' "$k" "${STEP_SCRIPT[$k]##*/}" "${STEP_FLAGS[$k]}"
    done
    cat <<EOF

Extra arguments are appended, e.g.:
  --step sample_qc --overwrite
  --step train_rf --model_id my_model

Config: set GNOMAD_QC_CONFIG to the path of your filled-in config.yaml,
and bind-mount that file into the container at that path.
EOF
}

if [[ "${1:-}" != "--step" || -z "${2:-}" ]]; then
    usage
    exit 1
fi

STEP="$2"
shift 2

SCRIPT="${STEP_SCRIPT[$STEP]:-}"
if [[ -z "$SCRIPT" ]]; then
    echo "Unknown step: '$STEP'" >&2
    usage
    exit 1
fi
FLAGS="${STEP_FLAGS[$STEP]:-}"

: "${GNOMAD_QC_CONFIG:?Set GNOMAD_QC_CONFIG to your config.yaml path and bind-mount it into the container}"

# FLAGS is intentionally unquoted so it word-splits into separate arguments.
exec python3 "/opt/gnomad_qc/${SCRIPT}" ${FLAGS} "$@"
