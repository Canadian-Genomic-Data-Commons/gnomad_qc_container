#!/usr/bin/env python3
import hail as hl
import atexit
import os
import re
import shutil
import sys
import pandas as pd
import numpy as np
from pathlib import Path
import traceback
import threading, time
import argparse
import fcntl

# === Configuration ===
BASE_DIR = "/projects/HostSeq/hostseq_smk/data"
STATS = "/projects/HostSeq/hostseq_jc/stats"
SCRATCH = "/projects/HostSeq/scratch/gnomad"
SAMTOOLS_HEADER = "/projects/HostSeq/hail/samtools_header"
JOINT_CALL_SAMPLES = "/projects/HostSeq/hostseq_jc/utils/lists_and_spreadsheets/joint-called-samples.txt"
BATCH_SIZE = 100
GVCF_BATCH_SIZE = 40
LOCAL_GVCF_BATCH_SIZE = 100
MIN_FILESIZE = 1_000_000
GATK_FOLD = "gatk-v4.2.5.0"
SUFFIX = ".g.vcf.gz"
RIGHT_LENGTH = 2934876451
LENGTH_SHEET = 15
NUM_ROWS = 25
CORR_RATIO = 32173469172 / RIGHT_LENGTH

# CLI arguments
parser = argparse.ArgumentParser()
parser.add_argument("--stragglers", action="store_true", help="Only process missing/unimported samples")
args = parser.parse_args()

# Determine array batch index early so we can use it for per-task local temp dirs
task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))
running_slurm = "SLURM_JOB_ID" in os.environ or "SLURM_ARRAY_TASK_ID" in os.environ
mode = "SLURM" if running_slurm else "LOCAL"
print(f"[task {task_id}] Running in {mode} mode")

# Prepare node-local temp for Spark (prefer SLURM_TMPDIR)
local_tmp = os.environ.get("SLURM_TMPDIR", "/tmp")
spark_tmp = os.path.join(local_tmp, f"spark_tmp_batch_{task_id}")
os.makedirs(spark_tmp, exist_ok=True)

# scratch dir for shared workspace and log
logpath = os.path.join(SCRATCH, f"combine_batch_{task_id}.log")
os.makedirs(SCRATCH, exist_ok=True)

# Make sure any running Hail is stopped before init
try:
    hl.stop()
except Exception:
    pass

print(f"[task {task_id}] spark_tmp={spark_tmp} logpath={logpath} scratch={SCRATCH}")

# Init hail with spark local dir pointed at node-local tmp (reduces NFS contention)
if mode == "SLURM":
    hl.init(
        log=logpath,
        spark_conf={
            # Resource settings
            "spark.driver.memory": "32g",
            "spark.executor.memory": "280g",
            "spark.executor.instances": "1",
            "spark.executor.cores": "48",
            "spark.task.cpus": "1",
            "spark.local.dir": spark_tmp,
            # Partitioning (more tasks, more parallelism)
            "spark.sql.files.maxPartitionBytes": "134217728",
            "spark.hadoop.mapreduce.input.fileinputformat.split.minsize": "134217728",
            # Stability
            "spark.executor.heartbeatInterval": "60s",
            "spark.network.timeout": "600s",
            # Adaptive query planning helps with skewed joins or large VCFs
            "spark.sql.adaptive.enabled": "true",
            "spark.sql.adaptive.skewJoin.enabled": "true",
        },
        default_reference="GRCh38",
    )
elif mode == "LOCAL":
    hl.init(
        log=logpath,
        spark_conf={
            # Resource settings
            "spark.driver.memory": "16g",
            "spark.executor.memory": "180g",
            "spark.executor.instances": "1",
            "spark.executor.cores": "32",
            "spark.task.cpus": "1",
            "spark.local.dir": spark_tmp,
            # Partitioning (more tasks, more parallelism)
            "spark.sql.files.maxPartitionBytes": "134217728",
            "spark.hadoop.mapreduce.input.fileinputformat.split.minsize": "134217728",
            # Stability
            "spark.executor.heartbeatInterval": "60s",
            "spark.network.timeout": "600s",
            # Adaptive query planning helps with skewed joins or large VCFs
            "spark.sql.adaptive.enabled": "true",
            "spark.sql.adaptive.skewJoin.enabled": "true",
        },
        default_reference="GRCh38",
    )


def log_heartbeat():
    while True:
        print(f"[task {task_id}] Heartbeat: still running at {time.ctime()}")
        time.sleep(300)


threading.Thread(target=log_heartbeat, daemon=True).start()

# Withdrawn / otherwise forbidden samples
FORBIDDEN_FILE = "/projects/HostSeq/hostseq_jc/utils/lists_and_spreadsheets/consent-withdrawn-clean.txt"
FORBIDDEN = set()
if os.path.exists(FORBIDDEN_FILE):
    with open(FORBIDDEN_FILE) as fh:
        FORBIDDEN = {line.strip() for line in fh if line.strip()}

DB_FILE = os.path.join(SCRATCH, "VDS_content.txt")
Path(DB_FILE).parent.mkdir(parents=True, exist_ok=True)
if not os.path.exists(DB_FILE):
    open(DB_FILE, "a").close()

ALREADY_IMPORTED = set()
ALREADY_IMPORTED_REFORMATTED = set()
if os.path.getsize(DB_FILE) > 0:
    try:
        df_db = pd.read_csv(DB_FILE, delimiter="\t", header=None, dtype=str)
        already_in = df_db[0].astype(str).tolist()
        for f in already_in:
            if "__" in f:
                f_trunc = f[11:].replace("_", "-")
                ALREADY_IMPORTED_REFORMATTED.add(f_trunc)
            else:
                ALREADY_IMPORTED_REFORMATTED.add(f.replace("_", "-"))
        ALREADY_IMPORTED.update(already_in)
    except Exception as e:
        print(f"[task {task_id}] Warning: Unable to read DB_FILE {DB_FILE}: {e}")
        ALREADY_IMPORTED = set()

all_samples = sorted(
    [s for s in os.listdir(BASE_DIR) if os.path.isdir(os.path.join(BASE_DIR, s))]
)
batch_start = task_id * BATCH_SIZE
batch_end = batch_start + BATCH_SIZE
if args.stragglers:
    df_jc = pd.read_csv(JOINT_CALL_SAMPLES, delimiter="\t", header=None, dtype=str) 
    jc_samples = df_jc[0].astype(str).tolist()
    stragglers = list(sorted(set(jc_samples) - set(already_in)))
    batch_samples = [
        s for s in stragglers
        if s not in FORBIDDEN and s not in ALREADY_IMPORTED and s not in ALREADY_IMPORTED_REFORMATTED
    ]
    batch_samples = stragglers[batch_start:batch_end]
    print(f"[task {task_id}] stragglers override mode: {len(batch_samples)} samples to process")
else:
    batch_start = task_id * BATCH_SIZE
    batch_end = batch_start + BATCH_SIZE
    batch_samples = all_samples[batch_start:batch_end]
    print(f"[task {task_id}] normal batching: samples {batch_start}:{batch_end} -> {len(batch_samples)} samples")


def human_bytes(n):
    # small helper for logs
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if n < 1024.0:
            return f"{n:.1f}{unit}"
        n /= 1024.0
    return f"{n:.1f}PiB"


def staged_gvcf_batches(gvcf_paths, batch_size, local_tmp, task_id, max_bytes=None):
    """
    Yield (local_batch_paths, local_batch_dir) for successive batches.

    Batches are built greedily, respecting both:
      - file-count limit `batch_size`
      - total bytes limit `max_bytes` (if provided)

    Each local_batch_dir is created under local_tmp and must be removed by caller.
    """
    batch_num = 0
    if max_bytes is None:
        # compute default: use 60% of free space on local_tmp with a floor of 20 GiB
        try:
            usage = shutil.disk_usage(local_tmp)
            free = usage.free
            calc = int(free * 0.6)
            floor = 20 * 1024**3
            max_bytes = max(calc, floor)
        except Exception:
            # fallback if disk_usage fails
            max_bytes = 20 * 1024**3
    print(
        f"[task {task_id}] LOCAL_GVCF_MAX_BYTES = {max_bytes} ({human_bytes(max_bytes)})"
    )

    i = 0
    n = len(gvcf_paths)
    while i < n:
        local_batch = []
        batch_bytes = 0
        local_batch_dir = os.path.join(local_tmp, f"gvcfs_batch_{task_id}_{batch_num}")
        os.makedirs(local_batch_dir, exist_ok=True)

        # Fill batch by file count and bytes
        while i < n and len(local_batch) < batch_size:
            src = gvcf_paths[i]
            try:
                sz = os.path.getsize(src)
            except Exception:
                sz = 0
            # If the file itself is larger than max_bytes and the batch is empty,
            # include it but warn. If batch not empty and this file would exceed
            # max_bytes, break to start a new batch.
            if len(local_batch) == 0:
                if sz > max_bytes:
                    print(
                        f"[task {task_id}] WARNING: single gVCF {src} size {human_bytes(sz)} exceeds batch byte limit {human_bytes(max_bytes)}; it will be processed alone."
                    )
                    local_batch.append(src)
                    batch_bytes += sz
                    i += 1
                    break
                else:
                    # fits into empty batch
                    local_batch.append(src)
                    batch_bytes += sz
                    i += 1
            else:
                if batch_bytes + sz > max_bytes:
                    # stop this batch and start a new one
                    break
                else:
                    local_batch.append(src)
                    batch_bytes += sz
                    i += 1

        print(
            f"[task {task_id}] preparing batch #{batch_num}: {len(local_batch)} files, total {human_bytes(batch_bytes)}"
        )
        local_batch_paths = []
        for gvcf in local_batch:
            local_path = os.path.join(local_batch_dir, os.path.basename(gvcf))
            if not os.path.exists(local_path):
                # copy to node-local area
                try:
                    shutil.copy2(gvcf, local_path)
                    print(f"[task {task_id}] copied {gvcf} -> {local_path}")
                except Exception as e:
                    print(f"[task {task_id}] ERROR copying {gvcf} -> {local_path}: {e}")
                    # clean up and re-raise so caller can inspect
                    try:
                        shutil.rmtree(local_batch_dir, ignore_errors=True)
                    except Exception:
                        pass
                    raise
            local_batch_paths.append(local_path)

        yield local_batch_paths, local_batch_dir
        batch_num += 1


def get_size(path):
    try:
        return os.path.getsize(path)
    except Exception:
        return 0


def passed_qc(sample):
    """
    Return True if sample passes all QC and should be included.
    NOTE: this function does NOT try to import the VCF (that was causing failures).
    """
    # Skip forbidden or already imported
    if sample in FORBIDDEN:
        print(f"[{sample}] skipped: in forbidden list")
        return False
    if sample in ALREADY_IMPORTED or sample in ALREADY_IMPORTED_REFORMATTED:
        print(f"[{sample}] skipped: already imported")
        return False

    # Check gVCF existence/size
    gvcf = os.path.join(BASE_DIR, sample, GATK_FOLD, sample + SUFFIX)
    if not os.path.exists(gvcf):
        print(f"[{sample}] skipped: gVCF missing ({gvcf})")
        return False
    if get_size(gvcf) < MIN_FILESIZE:
        print(f"[{sample}] skipped: gVCF too small ({get_size(gvcf)} bytes)")
        return False

    # Bypass if doing stragglers only
    if args.stragglers:
        return True

    # Stats file paths (centralized STATS dir)
    flagstat = os.path.join(STATS, "flagstats", sample + ".md.recal.bam.flagstat")
    bamstat = os.path.join(STATS, "bamstats", sample + ".md.recal.bam.bamstat")
    coverage = os.path.join(STATS, "coverage", sample + ".md.recal.bam.coverage")
    idxstat = os.path.join(STATS, "idxstats", sample + ".md.recal.bam.idxstat")
    input_file = os.path.join(STATS, "input", sample)

    for p, name in [
        (flagstat, "flagstat"),
        (bamstat, "bamstat"),
        (coverage, "coverage"),
        (idxstat, "idxstat"),
        (input_file, "input"),
    ]:
        if not os.path.exists(p):
            print(f"[{sample}] skipped: {name} missing ({p})")
            return False
        if get_size(p) == 0:
            print(f"[{sample}] skipped: {name} zero-size ({p})")
            return False

    # small, robust parsing checks (same as before)
    try:
        df_inp = pd.read_csv(input_file, delimiter="\t", header=None, nrows=NUM_ROWS)
        df_inp.iloc[:, 0] = (
            df_inp.iloc[:, 0]
            .astype(str)
            .str.replace(r"\D+", "", regex=True)
            .astype(int)
        )
    except Exception as e:
        print(f"[{sample}] skipped: error reading input file: {e}")
        return False

    try:
        df_is = pd.read_csv(idxstat, delimiter="\t", header=None, nrows=NUM_ROWS)
        if df_is.shape[1] < 4:
            print(f"[{sample}] skipped: idxstat unexpected format")
            return False
        df_is.columns = ["chr", "seqlen", "mappedseg", "unmappedseg"]
    except Exception as e:
        print(f"[{sample}] skipped: error reading idxstat: {e}")
        return False

    try:
        df_cov = pd.read_csv(coverage, delimiter="\t", nrows=NUM_ROWS)
    except Exception as e:
        print(f"[{sample}] skipped: error reading coverage: {e}")
        return False

    try:
        df_fs = pd.read_csv(
            flagstat, delimiter=" ", names=list(range(LENGTH_SHEET)), nrows=NUM_ROWS
        )
    except Exception as e:
        print(f"[{sample}] skipped: error reading flagstat: {e}")
        return False

    try:
        # check1: minimum count across chr1-22+X (first 23 entries)
        check1_val = np.min(df_inp.iloc[0:23, 0].values)
        check1 = check1_val < 250000

        # check2: mapped/unmapped ratio max across chr1-22
        mapped = df_is["mappedseg"].astype(int)
        unmapped = df_is["unmappedseg"].astype(int)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratios = np.where(mapped == 0, 0.0, unmapped / mapped)
        check2 = np.nanmax(ratios[0:21]) > 0.15

        # check3: coverage min across chr1-22
        if "coverage" in df_cov.columns:
            cov_col = "coverage"
        else:
            cov_col = df_cov.columns[0]
        covs = df_cov[cov_col].astype(float)
        check3 = np.nanmin(covs[0:21]) < 70

        # check4: total reads from flagstat first column
        if isinstance(df_fs.iloc[0, 0], str):
            counts = df_fs.iloc[:, 0].str.replace(r"\D+", "", regex=True).astype(int)
        else:
            counts = (
                pd.to_numeric(df_fs.iloc[:, 0], errors="coerce").fillna(0).astype(int)
            )
        check4 = counts.iloc[0] < 5.5e8

        # check5: X coverage from bamstat
        x_cov = None
        eff_size = ""
        with open(bamstat, "r") as fh:
            for line in fh:
                if line.startswith("Effective_Genome_size:"):
                    eff_size = line
                if line.startswith("Estimate_for_genome_X_coverage:"):
                    try:
                        s = line.strip()
                        m = re.search(r"[-+]?\d*\.?\d+([eE][-+]?\d+)?", s)
                        if m:
                            x_cov = float(m.group(0))
                    except Exception:
                        x_cov = None
        if x_cov is None:
            x_cov = 0.0
        if str(RIGHT_LENGTH) not in eff_size:
            x_cov = x_cov * CORR_RATIO
        check5 = x_cov < 10

        if any([check1, check2, check3, check4, check5]):
            print(
                f"[{sample}] failed QC checks: check1={check1_val}, check2_max={np.nanmax(ratios):.3f}, check3_min={np.nanmin(covs):.2f}, total_reads={counts.iloc[0]}, Xcov={x_cov:.2f}"
            )
            return False
    except Exception as e:
        print(f"[{sample}] skipped: error during QC calculations: {e}")
        return False
    return True


# Build list of gvcf paths that pass QC checks (no Hail import check)
gvcf_paths = []
for sample in batch_samples:
    gvcf = os.path.join(BASE_DIR, sample, GATK_FOLD, sample + SUFFIX)
    if passed_qc(sample):
        gvcf_paths.append(gvcf)

print(f"[task {task_id}] Batch {task_id}: {len(gvcf_paths)} gVCFs passed QC")

if not gvcf_paths:
    print(f"[task {task_id}] No valid gVCFs in this batch; exiting")
    sys.exit(0)

# Output paths (use .tmp first to avoid incomplete VDS confusion)
final_vds = os.path.join(SCRATCH, f"batch_{task_id}_stragglers.vds") if args.stragglers else os.path.join(SCRATCH, f"batch_{task_id}.vds")
temp_vds = final_vds + ".tmp"
temp_dir = os.path.join(SCRATCH, f"temp_batch_{task_id}")
os.makedirs(temp_dir, exist_ok=True)

# If a previous temp directories exist, clean up first (only under SCRATCH)
if os.path.exists(temp_vds):
    print(f"[task {task_id}] Found leftover temp VDS {temp_vds}; deleting")
    shutil.rmtree(temp_vds, ignore_errors=True)
if os.path.exists(temp_dir):
    print(f"[task {task_id}] Removing stale temp directory {temp_dir}")
    shutil.rmtree(temp_dir, ignore_errors=True)

os.makedirs(temp_dir, exist_ok=True)
print(f"[task {task_id}] Running hail.vds.new_combiner on {len(gvcf_paths)} gVCFs")

try:
    
    lock_file = final_vds + ".lock"
    with open(lock_file, "w") as lf:
        fcntl.flock(lf, fcntl.LOCK_EX)
        # Skip if final batch already exists
        if os.path.exists(final_vds):
            print(f"[task {task_id}] Batch VDS {final_vds} already exists; skipping combiner")
            sys.exit(0)

        combiner = hl.vds.new_combiner(
            output_path=temp_vds,
            temp_path=temp_dir,
            gvcf_paths=gvcf_paths,
            gvcf_sample_names=[
                os.path.basename(path).replace(SUFFIX, "") for path in gvcf_paths
            ],
            gvcf_external_header=SAMTOOLS_HEADER,
            gvcf_batch_size=GVCF_BATCH_SIZE if mode == "SLURM" else LOCAL_GVCF_BATCH_SIZE,
            reference_genome="GRCh38",
            use_genome_default_intervals=True,
        )
        combiner.run()

        # Move to final location
        os.rename(temp_vds, final_vds)
        fcntl.flock(lf, fcntl.LOCK_UN)
        print(f"[task {task_id}] Wrote final VDS to {final_vds}")

except Exception:
    print(f"[task {task_id}] FATAL: combiner.run() failed; temp_dir left at {temp_dir}")
    traceback.print_exc()
    try:
        hl.stop()
    except Exception as e:
        print(
            f"[task {task_id}] Warning: hl.stop() failed during exception handling: {e}"
        )
    try:
        if spark_tmp.startswith(local_tmp):
            shutil.rmtree(spark_tmp, ignore_errors=True)
            print(f"[task {task_id}] Removed spark_tmp {spark_tmp}")
    except Exception as e:
        print(f"[task {task_id}] Warning: could not remove spark_tmp {spark_tmp}: {e}")
    sys.exit(1)

try:
    hl.stop()
except Exception as e:
    print(f"[task {task_id}] Warning: hl.stop() failed: {e}")

try:
    if spark_tmp.startswith(local_tmp):
        shutil.rmtree(spark_tmp, ignore_errors=True)
        print(f"[task {task_id}] Removed node-local spark tmp {spark_tmp}")
except Exception as e:
    print(f"[task {task_id}] Warning: could not remove spark_tmp {spark_tmp}: {e}")


def cleanup():
    try:
        hl.stop()
    except:
        pass
    try:
        shutil.rmtree(spark_tmp, ignore_errors=True)
    except:
        pass


atexit.register(cleanup)

# Remove any leftover local gvcf batch directories under local_tmp (best-effort)
for maybe in os.listdir(local_tmp):
    if maybe.startswith(f"gvcfs_batch_{task_id}_"):
        path = os.path.join(local_tmp, maybe)
        try:
            shutil.rmtree(path, ignore_errors=True)
            print(f"[task {task_id}] Removed leftover local gvcf batch {path}")
        except Exception:
            pass

# Append samples to DB only after final VDS exists
try:
    with open(DB_FILE, "a+") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        for g in gvcf_paths:
            sample = os.path.basename(g).replace(SUFFIX, "")
            fh.write(sample + "\n")
        fcntl.flock(fh, fcntl.LOCK_UN)
    print(f"[task {task_id}] Appended {len(gvcf_paths)} samples to {DB_FILE}")
except Exception as e:
    print(f"[task {task_id}] Warning: unable to append to DB_FILE: {e}")

print(f"[task {task_id}] DONE")
