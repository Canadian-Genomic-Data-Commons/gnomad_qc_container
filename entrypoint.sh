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

# Ordered to follow the pipeline sequence in the gnomad_qc v3 README

STEP_ORDER=(
    bootstrap_meta
    create_last_END_positions
    sample_qc
    impute_sex
    compute_hard_filters
    generate_qc_annotations
    compute_qc_mt
    pc_relate
    run_pca
    assign_pops
    ancestry_classification
    calculate_inbreeding
    calculate_clinvar
    apply_regressed_filters
    compute_related_samples_to_drop
    find_dups
    generate_metadata
    apply_stratified_filters
    infer_families
    run_mendel_errors
    finalize_ped
    export_info_vcf
    load_vqsr
    freq_run_dense_annotations
    freq_correct_high_ab_hets
    generate_allele_data
    generate_ac
    generate_fam_stats
    annotate_for_rf
    train_rf
    apply_rf
    create_bin_ht
    create_aggregated_bin_ht
    extract_truth_samples
    merge_with_truth_data
    bin_truth_sample_concordance
    final_filter
    vep
    export_transmitted_singletons_vcf
    create_release_sites_ht
    create_hgdp_tgp_subset
    make_var_annot_hists
    prepare_vcf_data_release
    compute_coverage
    v2_pc_relate
)

declare -A STEP_SCRIPT=(
    [bootstrap_meta]="gnomad_qc/v3/resources/meta.py"
    [create_last_END_positions]="gnomad_qc/v3/load_data/create_last_END_positions.py"
    [sample_qc]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [impute_sex]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [compute_hard_filters]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [compute_qc_mt]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [pc_relate]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [run_pca]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [assign_pops]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [ancestry_classification]="gnomad_qc/v3/sample_qc/ancestry_classification.py"
    [calculate_inbreeding]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [calculate_clinvar]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [apply_stratified_filters]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [apply_regressed_filters]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [compute_related_samples_to_drop]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [generate_metadata]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [find_dups]="gnomad_qc/v3/sample_qc/create_fam.py"
    [infer_families]="gnomad_qc/v3/sample_qc/create_fam.py"
    [run_mendel_errors]="gnomad_qc/v3/sample_qc/create_fam.py"
    [finalize_ped]="gnomad_qc/v3/sample_qc/create_fam.py"
    [generate_qc_annotations]="gnomad_qc/v3/annotations/generate_qc_annotations.py"
    [export_info_vcf]="gnomad_qc/v3/annotations/generate_qc_annotations.py"
    [load_vqsr]="gnomad_qc/v3/load_data/load_vqsr.py"
    [annotate_for_rf]="gnomad_qc/v3/variant_qc/random_forest.py"
    [train_rf]="gnomad_qc/v3/variant_qc/random_forest.py"
    [apply_rf]="gnomad_qc/v3/variant_qc/random_forest.py"
    [create_bin_ht]="gnomad_qc/v3/variant_qc/evaluation.py"
    [create_aggregated_bin_ht]="gnomad_qc/v3/variant_qc/evaluation.py"
    [extract_truth_samples]="gnomad_qc/v3/variant_qc/evaluation.py"
    [merge_with_truth_data]="gnomad_qc/v3/variant_qc/evaluation.py"
    [bin_truth_sample_concordance]="gnomad_qc/v3/variant_qc/evaluation.py"
    [final_filter]="gnomad_qc/v3/variant_qc/final_filter.py"
    [freq_run_dense_annotations]="gnomad_qc/v4/annotations/generate_freq.py"
    [freq_correct_high_ab_hets]="gnomad_qc/v4/annotations/generate_freq.py"
    [generate_allele_data]="gnomad_qc/v3/annotations/generate_qc_annotations.py"
    [generate_ac]="gnomad_qc/v3/annotations/generate_qc_annotations.py"
    [generate_fam_stats]="gnomad_qc/v3/annotations/generate_qc_annotations.py"
    [vep]="gnomad_qc/v3/annotations/generate_qc_annotations.py"
    [export_transmitted_singletons_vcf]="gnomad_qc/v3/annotations/generate_qc_annotations.py"
    [create_release_sites_ht]="gnomad_qc/v3/create_release/create_release_sites_ht.py"
    [create_hgdp_tgp_subset]="gnomad_qc/v3/create_release/create_hgdp_tgp_subset.py"
    [make_var_annot_hists]="gnomad_qc/v3/create_release/make_var_annot_hists.py"
    [prepare_vcf_data_release]="gnomad_qc/v3/create_release/prepare_vcf_data_release.py"
    [compute_coverage]="gnomad_qc/v3/load_data/compute_coverage.py"
    [v2_pc_relate]="gnomad_qc/v3/sample_qc/v2_pc_relate.py"
)

# Empty means the script takes no step-selecting flags (it runs top to bottom).
declare -A STEP_FLAGS=(
    [bootstrap_meta]=""
    [create_last_END_positions]=""
    [sample_qc]="--sample_qc"
    [impute_sex]="--impute_sex"
    [compute_hard_filters]="--compute_hard_filters"
    [compute_qc_mt]="--compute_qc_mt"
    [pc_relate]="--run_pc_relate"
    [run_pca]="--run_pca"
    [assign_pops]="--assign_pops"
    [ancestry_classification]=""
    [calculate_inbreeding]="--calculate_inbreeding"
    [calculate_clinvar]="--calculate_clinvar"
    [apply_stratified_filters]="--apply_stratified_filters"
    [apply_regressed_filters]="--apply_regressed_filters"
    [compute_related_samples_to_drop]="--compute_related_samples_to_drop"
    [generate_metadata]="--generate_metadata"
    [find_dups]="--find_dups"
    [infer_families]="--infer_families"
    [run_mendel_errors]="--run_mendel_errors"
    [finalize_ped]="--finalize_ped"
    [generate_qc_annotations]="--compute_info --split_info"
    [export_info_vcf]="--export_info_vcf"
    [load_vqsr]=""
    [annotate_for_rf]="--annotate_for_rf --impute_features"
    [train_rf]="--train_rf"
    [apply_rf]="--apply_rf"
    [create_bin_ht]="--create_bin_ht"
    [create_aggregated_bin_ht]="--create_aggregated_bin_ht"
    [extract_truth_samples]="--extract_truth_samples"
    [merge_with_truth_data]="--merge_with_truth_data"
    [bin_truth_sample_concordance]="--bin_truth_sample_concordance"
    [final_filter]=""
    [freq_run_dense_annotations]="--run-freq-and-dense-annotations"
    [freq_correct_high_ab_hets]="--correct-for-high-ab-hets"
    [generate_allele_data]="--generate_allele_data"
    [generate_ac]="--generate_ac"
    [generate_fam_stats]="--generate_fam_stats"
    [vep]="--vep"
    [export_transmitted_singletons_vcf]="--export_transmitted_singletons_vcf"
    [create_release_sites_ht]=""
    [create_hgdp_tgp_subset]=""
    [make_var_annot_hists]=""
    [prepare_vcf_data_release]=""
    [compute_coverage]=""
    [v2_pc_relate]="--run_pc_relate"
)

# Steps that need arguments you must supply yourself (which model, which
# cutoffs). Without them the script exits with an argument error.
declare -A STEP_NEEDS=(
    [apply_rf]="--model_id <id>        e.g. --model_id rf_c4e89eb8"
    [create_bin_ht]="--model_id <id>        e.g. --model_id rf_c4e89eb8"
    [create_aggregated_bin_ht]="--model_id <id>        e.g. --model_id rf_c4e89eb8"
    [extract_truth_samples]="--model_id <id>"
    [merge_with_truth_data]="--model_id <id>"
    [bin_truth_sample_concordance]="--model_id <id>"
    [load_vqsr]="--vqsr_vcf_path <path>   plus --header_path for the fixed VQSR header"
    [final_filter]="--model_id and cutoffs  e.g. --model_id rf_c4e89eb8 --model_name RF --score_name RF --snp_bin_cutoff 90 --indel_bin_cutoff 80"
)

# Steps that run on the host rather than in this container. VEP runs in its own
# container, which can't be started from inside this one; host/run_vep.sh runs
# the step with the same code and config.yaml.
declare -A STEP_RUN_ON_HOST=(
    [vep]="host/run_vep.sh"
)

# Steps that can't run outside the Broad: they read or write gnomAD's private
# storage (gs://gnomad/ and its temp buckets). To use one, change its resource
# paths to read from config.yaml, as the other steps do.
declare -A STEP_UNSUPPORTED=(
    [calculate_inbreeding]="reads qc.mt() from gs://gnomad/"
    [extract_truth_samples]="checkpoints to gs://gnomad-tmp-4day/"
    [merge_with_truth_data]="truth data under gs://gnomad/"
    [bin_truth_sample_concordance]="truth data under gs://gnomad/"
    [export_transmitted_singletons_vcf]="writes under gs://gnomad/"
    [make_var_annot_hists]="release data under gs://gnomad/"
    [prepare_vcf_data_release]="release data under gs://gnomad/"
    [v2_pc_relate]="reads gs://gnomad-tmp/v2_exomes_v3_joint_qc.mt"
)

usage() {
    cat <<EOF
Usage: --step <name> [extra args passed through to the script]

Steps, in SOP pipeline order (VDS creation happens before these, via build_vds/):

EOF
    for k in "${STEP_ORDER[@]}"; do
        printf '  %-28s %s %s\n' "$k" "${STEP_SCRIPT[$k]##*/}" "${STEP_FLAGS[$k]}"
        if [[ -n "${STEP_UNSUPPORTED[$k]:-}" ]]; then
            printf '  %-28s   ^ UNSUPPORTED here: %s\n' "" "${STEP_UNSUPPORTED[$k]}"
        fi
        if [[ -n "${STEP_RUN_ON_HOST[$k]:-}" ]]; then
            printf '  %-28s   ^ runs on the host, not in this container: %s\n' "" "${STEP_RUN_ON_HOST[$k]}"
        fi
        if [[ -n "${STEP_NEEDS[$k]:-}" ]]; then
            printf '  %-28s   ^ also needs: %s\n' "" "${STEP_NEEDS[$k]}"
        fi
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
if [[ -n "${STEP_RUN_ON_HOST[$STEP]:-}" ]]; then
    cat >&2 <<MSG
Step '$STEP' runs on the host, not inside this container.
VEP runs in its own container, and one container cannot start another here.
Run it with ${STEP_RUN_ON_HOST[$STEP]} from the gnomad_qc_container repository,
using the same config.yaml; it writes vep.ht to write_dir for later steps.
MSG
    exit 1
fi

if [[ -n "${STEP_UNSUPPORTED[$STEP]:-}" ]]; then
    cat >&2 <<MSG
Step '$STEP' is not supported outside the Broad: ${STEP_UNSUPPORTED[$STEP]}.
That storage is not public, so a local resource mirror cannot substitute for it.
To enable this step, convert its resource paths to read from config.yaml, as the
rest of the pipeline does, then remove it from STEP_UNSUPPORTED in this script.
MSG
    exit 1
fi

FLAGS="${STEP_FLAGS[$STEP]:-}"

: "${GNOMAD_QC_CONFIG:?Set GNOMAD_QC_CONFIG to your config.yaml path and bind-mount it into the container}"

exec python3 "/opt/gnomad_qc/${SCRIPT}" ${FLAGS} "$@"
