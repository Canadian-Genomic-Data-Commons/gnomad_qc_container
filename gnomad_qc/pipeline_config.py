"""Shared Spark/Hail configuration for gnomad_qc pipeline scripts.

Every pipeline step reads the same config.yaml (path set via the
GNOMAD_QC_CONFIG environment variable) instead of hardcoding its own paths
and Spark settings, so a lab can point the container at their own storage
and cluster sizing without editing pipeline code.
"""
import os
import sys

import hail as hl
import yaml

DEFAULT_CONFIG_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "v3",
    "configuration",
    "config.yaml",
)
CONFIG_PATH = os.environ.get("GNOMAD_QC_CONFIG", DEFAULT_CONFIG_PATH)

with open(CONFIG_PATH, "r") as f:
    _CFG = yaml.load(f, Loader=yaml.FullLoader) or {}

_PATHS = _CFG.get("paths") or {}
_SPARK = _CFG.get("spark") or {}
_JAVA = _CFG.get("java") or {}
_STEP_OVERRIDES = _CFG.get("step_overrides") or {}
_REFERENCE = _CFG.get("reference") or {}

BASE = _PATHS.get("base_dir")
WRITE_DIR = _PATHS.get("write_dir")
LOG_DIR = _PATHS.get("log_dir")
CODE_DIR = _PATHS.get("code_dir")
VDS_FILE = _PATHS.get("vds_file")

GRCH38_FASTA = _REFERENCE.get("grch38_fasta")
GRCH38_FASTA_INDEX = _REFERENCE.get("grch38_fasta_index")
GRCH37_FASTA = _REFERENCE.get("grch37_fasta")
GRCH37_FASTA_INDEX = _REFERENCE.get("grch37_fasta_index")

if not BASE:
    raise ValueError(
        f"'paths.base_dir' is not set in {CONFIG_PATH}. Point GNOMAD_QC_CONFIG at "
        "your filled-in config.yaml before running a pipeline step."
    )

# Only needed outside the container, where gnomad_qc isn't already on
# PYTHONPATH. The Dockerfile sets PYTHONPATH itself, so this is a no-op there.
if CODE_DIR:
    sys.path.insert(0, os.path.join(CODE_DIR, "gnomad_qc"))

if _JAVA.get("tool_options") and "JAVA_TOOL_OPTIONS" not in os.environ:
    os.environ["JAVA_TOOL_OPTIONS"] = _JAVA["tool_options"]

_SPARK_KEY_MAP = {
    "driver_memory": "spark.driver.memory",
    "executor_memory": "spark.executor.memory",
    "executor_instances": "spark.executor.instances",
    "executor_cores": "spark.executor.cores",
    "dynamic_allocation": "spark.dynamicAllocation.enabled",
    "shuffle_service": "spark.shuffle.service.enabled",
    "local_dir": "spark.local.dir",
    "memory_fraction": "spark.memory.fraction",
}


def _spark_conf(step_name: str) -> dict:
    """Merge global `spark` settings with `step_overrides[step_name]`."""
    settings = dict(_SPARK)
    settings.update(_STEP_OVERRIDES.get(step_name) or {})
    conf = {}
    for key, spark_key in _SPARK_KEY_MAP.items():
        value = settings.get(key)
        if value not in (None, ""):
            conf[spark_key] = str(value)
    return conf


def init_hail(step_name: str, default_reference: str = "GRCh38") -> None:
    """
    Stop any existing Hail context and start one using the shared config.

    :param step_name: Name of the pipeline step, used for the log file name
        and to look up `step_overrides` in config.yaml (e.g.
        "ancestry_classification").
    :param default_reference: Reference genome passed to hl.init.
    :return: None.
    """
    hl.stop()
    log_path = os.path.join(BASE, LOG_DIR, f"{step_name}.log") if LOG_DIR else None
    hl.init(
        log=log_path,
        spark_conf=_spark_conf(step_name),
        default_reference=default_reference,
    )
