#!/usr/bin/env python3
import os
import glob
import math
import subprocess

SCRATCH_DIR = "/projects/HostSeq/scratch/gnomad"
FINAL_VDS = "/projects/HostSeq/hail/hostseq.vds"
MERGE_FANIN = 10


def submit_merge_jobs(level, num_inputs):
    num_jobs = math.ceil(num_inputs / MERGE_FANIN)
    all_done = True
    for i in range(num_jobs):
        output_vds = f"{SCRATCH_DIR}/merge_{level}_{i}.vds"
        if not os.path.exists(output_vds):
            all_done = False
            break
    if all_done:
        print(f":white_check_mark: Level {level} merges already exist, skipping")
        return num_jobs

    print(
        f":rocket: Submitting level {level} merge ({num_inputs} inputs → {num_jobs} jobs)"
    )
    cmd = [
        "sbatch",
        f"--array=0-{num_jobs-1}",
        f"--export=MERGE_LEVEL={level}",
        "merge_batches.sh",
    ]
    subprocess.run(cmd, check=True)
    return num_jobs


def main():
    batch_inputs = glob.glob(f"{SCRATCH_DIR}/batch_*.vds")
    if not batch_inputs:
        raise RuntimeError("No batch VDS files found. Run combine_batch.py first.")

    num_inputs = len(batch_inputs)
    level = 1

    while num_inputs > 1:
        num_outputs = submit_merge_jobs(level, num_inputs)
        num_inputs = num_outputs
        level += 1

    print(f":tada: Final merge will write to {FINAL_VDS}")


if __name__ == "__main__":
    main()
