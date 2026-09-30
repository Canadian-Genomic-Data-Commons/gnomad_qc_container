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

# gnomad_methods resolves its public resources to gs:// URLs by default. On a
# cluster whose compute nodes have no internet, that makes Hail's GCS connector
# query the cloud metadata service at 169.254.169.254 and hang until it times
# out. Pointing the library at a local mirror avoids the network entirely: it
# rewrites gs://<public-bucket>/<rest> to <dir>/<rest>, so the mirror must keep
# the directory structure that follows the bucket name.
PUBLIC_RESOURCE_DIR = _REFERENCE.get("public_resource_dir")
if PUBLIC_RESOURCE_DIR and "GNOMAD_DEFAULT_PUBLIC_RESOURCE_SOURCE" not in os.environ:
    os.environ["GNOMAD_DEFAULT_PUBLIC_RESOURCE_SOURCE"] = PUBLIC_RESOURCE_DIR

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


def _hail_already_running() -> bool:
    """Whether a Hail context exists already (so the JVM is up)."""
    try:
        from hail.utils.java import Env

        return Env._hc is not None
    except Exception:
        return False


def init_hail(step_name: str, default_reference: str = "GRCh38") -> None:
    """
    Start Hail using the shared config.

    :param step_name: Name of the pipeline step, used for the log file name
        and to look up `step_overrides` in config.yaml (e.g.
        "ancestry_classification").
    :param default_reference: Reference genome passed to hl.init.
    :return: None.
    :raises RuntimeError: If Hail is already running, since the configured
        Spark settings could not take effect.
    """
    # A JVM's heap is fixed at launch, so hl.stop() + hl.init() cannot raise
    # driver memory. If something started Hail first (importing a module that
    # touches Hail at import time, say), the configured memory is silently
    # discarded -- a step asking for 200g would quietly run in ~1g and die far
    # into a long job. Fail immediately with an explanation instead.
    if _hail_already_running():
        raise RuntimeError(
            f"Hail was already initialized before init_hail('{step_name}') ran, so "
            "the Spark settings in your config.yaml cannot be applied (a JVM's "
            "heap size is fixed when it launches). Call init_hail() before "
            "importing anything that touches Hail at import time, such as the "
            "gnomad_qc.v3/v4.resources modules."
        )
    hl.stop()
    log_path = os.path.join(BASE, LOG_DIR, f"{step_name}.log") if LOG_DIR else None
    hl.init(
        log=log_path,
        spark_conf=_spark_conf(step_name),
        default_reference=default_reference,
    )
