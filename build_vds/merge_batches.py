#!/usr/bin/env python3
import hail as hl
import argparse
import os
import shutil

# Shared scratch for logs and workspace (adjust as needed)
SCRATCH = "/projects/HostSeq/hostseq_jc/scratch"

def safe_rmtree(path, scratch_root=SCRATCH):
    """
    Safely remove a directory only if it's inside the scratch area.
    Prevents accidental deletion of global /tmp or other important paths.
    """
    abs_path = os.path.abspath(path)
    if abs_path.startswith(os.path.abspath(scratch_root)):
        print(f"🧹 Removing directory {abs_path}")
        shutil.rmtree(abs_path, ignore_errors=True)
    else:
        print(f"⚠️ Skipping deletion of {abs_path} (not inside scratch root)")

def main():
    parser = argparse.ArgumentParser(
        description="Merge multiple batch VDSs into a single VDS."
    )
    parser.add_argument(
        "--input-list", "-i", required=True, help="Text file with batch VDS paths (one per line)"
    )
    parser.add_argument("--output", "-o", required=True, help="Output VDS path")
    parser.add_argument(
        "--temp-dir", "-t", default=None,
        help="Optional temporary directory for merging (auto-created per task if not provided)"
    )
    args = parser.parse_args()

    if os.path.exists(args.output):
        print(f"✅ Merge output {args.output} already exists, skipping")
        return

    # Determine array batch index early
    task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))

    # Node-local Spark temp
    local_tmp = os.environ.get("SLURM_TMPDIR", "/tmp")
    spark_tmp = os.path.join(local_tmp, f"spark_tmp_merge_{task_id}")
    os.makedirs(spark_tmp, exist_ok=True)

    # Log path
    logpath = os.path.join(SCRATCH, f"merge_batches_{task_id}.log")

    # Use temp_vds to avoid partial output confusion
    temp_vds = args.output + ".tmp"

    # If temp_dir not provided, auto-create a per-task scratch subdir
    if args.temp_dir is None:
        args.temp_dir = os.path.join(SCRATCH, f"temp_merge_{task_id}")
    os.makedirs(args.temp_dir, exist_ok=True)

    # Stop any running Hail
    try:
        hl.stop()
    except Exception:
        pass

    # Init Hail
    hl.init(
        log=logpath,
        spark_conf={
            "spark.driver.memory": "64g",
            "spark.executor.memory": "64g",
            "spark.local.dir": spark_tmp,
            "spark.sql.files.maxPartitionBytes": "268435456",
            "spark.hadoop.mapreduce.input.fileinputformat.split.minsize": "268435456",
        },
        default_reference="GRCh38",
        tmp_dir=args.temp_dir,
    )

    # Read input list
    with open(args.input_list, "r") as f:
        batch_vds_paths = [line.strip() for line in f if line.strip()]

    if not batch_vds_paths:
        print("❌ No batch VDS paths found in input list")
        return

    # Clean leftovers safely
    if os.path.exists(temp_vds):
        print(f"[task {task_id}] Found leftover temp VDS {temp_vds}; deleting")
        safe_rmtree(temp_vds)
    safe_rmtree(args.temp_dir)

    # Combine VDS shards
    combiner = hl.vds.new_combiner(
        output_path=temp_vds,
        temp_path=args.temp_dir,
        vds_paths=batch_vds_paths,
        reference_genome="GRCh38",
        gvcf_batch_size=len(batch_vds_paths),  # combine all shards at once
        use_genome_default_intervals=True,
    )
    combiner.run()

    # Rename temp VDS to final only on success
    os.rename(temp_vds, args.output)
    print(f"✅ Merged {len(batch_vds_paths)} batches into {args.output}")

if __name__ == "__main__":
    main()