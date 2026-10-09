#!/bin/bash
# Run the `vep` step on the host, outside the gnomad_qc container.
#
# VEP runs in its own container (e.g. vep114.sif), and starting one container
# from inside another is generally not possible on HPC, so this step can't be
# dispatched like the others. It uses the same pipeline code and the same
# config.yaml, and writes vep.ht to write_dir, where create_release_sites_ht
# (run in the container as usual) picks it up.
#
# Needs on the host:
#   - a Python environment with Hail 0.2.134, gnomad (gnomad_methods), pyyaml
#     and slackclient -- the same stack as the image
#   - config.yaml with reference.vep_config pointing at a Hail VEP config JSON,
#     whose "command" launches VEP (e.g. a wrapper running your VEP container)
#
# Usage:
#   sbatch host/run_vep.sh /path/to/config.yaml /path/to/python/env
#
# Resources: VEP on a full genome callset took ~24h on 40 cores / 75G in the
# production benchmark. Adjust the #SBATCH lines to your callset.

#SBATCH --job-name=gnomad_qc_vep
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=150G
#SBATCH --time=48:00:00

set -euo pipefail

CONFIG=${1:?usage: run_vep.sh /path/to/config.yaml /path/to/python/env}
PYENV=${2:?usage: run_vep.sh /path/to/config.yaml /path/to/python/env}

# Repository root, so `import gnomad_qc.pipeline_config` resolves to this
# checkout. Under sbatch the script is copied to a spool directory, so prefer
# the directory it was submitted from.
REPO=${GNOMAD_QC_REPO:-$(cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}" && pwd)}
if [ ! -f "${REPO}/gnomad_qc/pipeline_config.py" ]; then
    echo "Can't find the repository at ${REPO}. Submit from the repository root," \
         "or set GNOMAD_QC_REPO." >&2
    exit 1
fi

module load python java gcc 2>/dev/null || true
source "${PYENV}/bin/activate"

export GNOMAD_QC_CONFIG="${CONFIG}"
export PYTHONPATH="${REPO}${PYTHONPATH:+:${PYTHONPATH}}"
export JAVA_TOOL_OPTIONS="-XX:+UseG1GC -XX:G1HeapRegionSize=32M"
# Node-local scratch for Spark and Hail, if the config uses ${SLURM_TMPDIR}.
[ -n "${SLURM_TMPDIR:-}" ] && mkdir -p "${SLURM_TMPDIR}/spark" "${SLURM_TMPDIR}/hail"

python "${REPO}/gnomad_qc/v3/annotations/generate_qc_annotations.py" --vep --overwrite
