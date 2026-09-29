# gnomad_qc_container

Containerized [gnomAD QC pipeline](https://github.com/broadinstitute/gnomad_qc) for
running on HPC clusters (Slurm + Apptainer/Singularity) instead of Google Dataproc.

Derived from the Broad Institute's `gnomad_qc`, with the `gnomad_methods` library
included as a pinned submodule. Site-specific paths and Spark settings have been
factored out of the pipeline scripts into a single `config.yaml`, so each site
points the container at its own storage and cluster sizing without editing code.

## Setup

Clone with submodules -- `gnomad_methods` is required to build the image:

```
git clone --recursive <repo-url>
cd gnomad_qc_container
```

If you already cloned without `--recursive`:

```
git submodule update --init gnomad_methods
```

## Configure

Copy the template and fill in your own values:

```
cp gnomad_qc/v3/configuration/config.yaml my_config.yaml
```

| Section | What to set |
| --- | --- |
| `paths` | `base_dir`, `write_dir`, `log_dir`, `vds_file` -- absolute paths to your data and output locations |
| `spark` | Driver/executor memory and cores, matched to the resources your scheduler allocates the job |
| `java` | JVM options; the G1GC defaults suit large-heap Hail workloads |
| `reference` | GRCh38 FASTA + `.fai` paths (no public default -- use your own copy) |
| `step_overrides` | Optional per-step Spark overrides, since steps range from ~10G to ~300G peak memory |

## Build

Requires Docker, so build on a workstation or VM -- HPC login nodes do not have it.

```
docker build -t gnomad_qc:latest .
```

Then convert for the cluster:

```
docker save gnomad_qc:latest -o gnomad_qc.tar
# transfer gnomad_qc.tar to the cluster, then:
module load apptainer
apptainer build gnomad_qc.sif docker-archive://gnomad_qc.tar
```

## Run

Each pipeline step is selected with `--step`. Bind-mount your data and config, and
point `GNOMAD_QC_CONFIG` at the config's path *inside* the container:

```
apptainer run \
    --bind /path/to/data:/data \
    --bind /path/to/my_config.yaml:/config/config.yaml \
    --env GNOMAD_QC_CONFIG=/config/config.yaml \
    gnomad_qc.sif --step ancestry_classification
```

Run without arguments to list the available steps:

```
apptainer run gnomad_qc.sif
```

Hail initialization is shared: every step calls `init_hail()` from
`gnomad_qc/pipeline_config.py`, which applies your Spark settings and writes a
per-step log into `log_dir`.

## Resource requirements

Benchmarked on ~4400 samples on Narval/DRAC. Steps vary enormously, so size each
Slurm job individually and mirror the values in `step_overrides`.

| Step | Cores | Peak memory | Est. time (h) |
| --- | --- | --- | --- |
| *(VDS creation -- see `build_vds/`, run before these)* | 40 | 300G | 120 |
| `create_last_END_positions` | 40 | 50G | 24 |
| `generate_qc_annotations` | 40 | 150G | 36 |
| `impute_sex` | 40 | 250G | 24 |
| `sample_qc` | 40 | 50G | 12 |
| `pc_relate` | 40 | 125G | 4 |
| `ancestry_classification` | 40 | 250G | 6 |
| `generate_metadata` | 40 | 10G | 0.5 |
| `find_dups` / `infer_families` / `run_mendel_errors` / `finalize_ped` | 40 | 100G | 4 |
| `annotate_for_rf` / `train_rf` / `apply_rf` | 64 | 200G | 6 |
| `create_bin_ht` / `create_aggregated_bin_ht` | 40 | 110G | 1 |
| `final_filter` | 40 | 75G | 1 |
| `freq_*` (split VDS, dense annotations, high-AB hets, finalize) | 64 | 250G | 15 |
| `vep` | 40 | 75G | 24 |
| `create_release_sites_ht` | 40 | 75G | 0.5 |

Total pipeline runtime is roughly 7-8 days. Note that several steps were
substantially over-allocated during benchmarking (`sample_qc` used 16.8G of 400G;
`create_last_END_positions` used 11.6G of 92G), so the memory column is an upper
bound, not a target.

Other guidance from the benchmarking work:

- **Single node is preferred.** Multi-node showed poor parallelization -- a 3-node
  `compute_qc_mt` run managed only 29% CPU efficiency.
- **Use fast local storage for Spark temp dirs** (`spark.local_dir`). Lustre is slow
  at the many-small-file access patterns Hail generates.
- **Do not set `PYSPARK_SUBMIT_ARGS`** -- it overrides `hl.init()` and causes memory
  conflicts. This container drives Spark through `init_hail()` instead.
- **Resuming:** steps accept `--overwrite`; scripts that checkpoint can resume from
  existing intermediates rather than restarting.

## Storage

| Category | Size |
| --- | --- |
| VDS files | 5.2TB |
| Intermediate checkpoints and scratch | 10-15TB |
| QC outputs (HTs, TSVs, metrics) | ~50GB |
| Long-term archive | 4-10TB |
| **Total recommended** | **20-30TB** (minimum 10TB) |

Peak usage runs 3-4x input data size during processing.
