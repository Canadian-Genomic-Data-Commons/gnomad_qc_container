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

Steps vary widely -- roughly 10G to 300G peak memory and 40-64 cores, with full
pipeline runtimes measured in days on a ~4400-sample callset. Size each Slurm job
per step and set `step_overrides` to match.
