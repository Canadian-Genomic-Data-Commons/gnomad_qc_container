#!/usr/bin/env bash
# Entrypoint for the gnomad_qc container. Maps a --step name to the
# corresponding pipeline script and runs it.
#
#   apptainer run --bind /path/to/data:/data gnomad_qc.sif \
#       --step ancestry_classification
#
# GNOMAD_QC_CONFIG must point to a filled-in config.yaml (paths + Spark
# settings) that is bind-mounted into the container -- see
# gnomad_qc/v3/configuration/config.yaml for the template/schema.
set -euo pipefail

declare -A STEPS=(
    [sample_qc]="gnomad_qc/v3/sample_qc/sample_qc.py"
    [ancestry_classification]="gnomad_qc/v3/sample_qc/ancestry_classification.py"
    [create_fam]="gnomad_qc/v3/sample_qc/create_fam.py"
    [make_pca_plot]="gnomad_qc/v3/sample_qc/make_pca_plot.py"
    [subpop_analysis]="gnomad_qc/v3/sample_qc/subpop_analysis.py"
    [v2_pc_relate]="gnomad_qc/v3/sample_qc/v2_pc_relate.py"
)

usage() {
    cat <<EOF
Usage: --step <name> [extra args passed through to the script]

Available steps:
$(for k in "${!STEPS[@]}"; do echo "  $k"; done | sort)

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

SCRIPT="${STEPS[$STEP]:-}"
if [[ -z "$SCRIPT" ]]; then
    echo "Unknown step: '$STEP'" >&2
    usage
    exit 1
fi

: "${GNOMAD_QC_CONFIG:?Set GNOMAD_QC_CONFIG to your config.yaml path and bind-mount it into the container}"

exec python3 "/opt/gnomad_qc/${SCRIPT}" "$@"
