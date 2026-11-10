# gnomad_qc_container

Containerized [gnomAD QC pipeline](https://github.com/broadinstitute/gnomad_qc) for
running on HPC clusters (Slurm + Apptainer/Singularity).

Derived from the Broad Institute's `gnomad_qc`, with the `gnomad_methods` library
included as a pinned submodule. Site-specific paths and Spark settings have been
factored out of the pipeline scripts into a single `config.yaml`, so each site
points the container at its own storage and cluster sizing without editing code.

The pipeline starts from a Hail VDS of your callset (GRCh38) and ends with a
release sites table, `release_sites.ht`.

## Setup

Clone the repository for the config template and the `host/` scripts. Use
`--recursive` if you plan to build the image yourself, since the build needs the
`gnomad_methods` submodule:

```
git clone --recursive https://github.com/Canadian-Genomic-Data-Commons/gnomad_qc_container.git
cd gnomad_qc_container
```

If you already cloned without `--recursive`:

```
git submodule update --init gnomad_methods
```

## Inputs

| Input | Where it goes |
| --- | --- |
| Hail VDS of your callset (GRCh38), written with Hail 0.2.134 or earlier (file format 1.7.0 or lower) | `paths.vds_file` in the config |
| GRCh38 FASTA and `.fai` | `reference.grch38_fasta` / `grch38_fasta_index` |
| Sample metadata TSV, named `meta_test.tsv` | Your `write_dir`, before the first step |
| VEP install and cache | Used by the `vep` step on the host (see [VEP](#vep)) |
| gnomAD reference tables | `base_dir/hail/`, `base_dir` and `reference.public_resource_dir` (see [Reference tables](#reference-tables)) |

The container runs Hail 0.2.134, which can read files in older formats but not
newer ones, so a VDS written with a later Hail version may not load. To match
exactly, build it with `pip install hail==0.2.134` (requires Java 11).

`meta_test.tsv` is tab-separated with a header row and one row per sample. The
`bootstrap_meta` step turns it into `meta.ht`. Required columns:

| Column | Meaning |
| --- | --- |
| `s` | Sample ID, exactly as in the VDS |
| `exclude` | `True` to drop the sample from QC |
| `releasable` | `True` if the sample can be included in the release |

Any other columns (e.g. `age`, a project population label) are carried into the
sample metadata as `project_meta` fields.

An existing `meta.ht` is never overwritten. After editing `meta_test.tsv`,
delete `meta.ht` and rerun `bootstrap_meta`.

### Reference tables

The pipeline also reads gnomAD reference tables: centromere and telomere
regions, truth sites, ancestry models, dbSNP and others. All are in gnomAD's
free public bucket (`gs://gcp-public-data--gnomad`); no Google account or
requester-pays billing is needed. `host/download_resources.sh` downloads them
with `gsutil` (part of the Google Cloud SDK). Run it on a node with internet
access:

```
nohup bash host/download_resources.sh /path/to/base_dir /path/to/gnomad_public_mirror > download_resources.log 2>&1 &
```

If it is interrupted, rerun it; only missing files are downloaded. When it
finishes, set `reference.public_resource_dir` to the mirror folder.

| Table | Where | Used by | Size |
| --- | --- | --- | --- |
| `gnomad.genomes.r3.0.sites.ht` | `base_dir` | `compute_qc_mt` (choosing QC sites) | ~340 GB |
| `purcell5k.ht` | `base_dir` | `compute_qc_mt` | <1 MB |
| `hg38.telomeresAndMergedCentromeres.ht` | `base_dir/hail` | sample QC, random forest, evaluation, `final_filter` | 2 MB |
| `LCRFromHengHg38.ht`, `GRCh38_segdups.ht` | `base_dir/hail` | `compute_qc_mt`, `create_release_sites_ht` | 35 MB |
| `truth_sites.ht`, `clinvar_20190923.ht` | `base_dir/hail` | random forest, evaluation | 530 MB |
| `gnomad.v3.1.pca_loadings.ht`, `gnomad.v3.1.RF_fit.onnx` | `base_dir/hail` | `ancestry_classification` | 490 MB |
| `dbsnp_b151_grch38_all_20180418.ht` | `base_dir/hail` | `create_release_sites_ht` | 19 GB |

The pipeline scripts read these from `base_dir` and `base_dir/hail`. The
`gnomad` library looks for some of the same tables under
`reference.public_resource_dir`, at the same paths as in gnomAD's bucket (e.g.
`resources/grch38/clinvar/clinvar_20190923.ht`); without them it tries gnomAD's
requester-pays bucket and fails. The script links those mirror paths to the
copies in `base_dir`, so nothing is downloaded twice. Bind the mirror folder
into the container as well as `base_dir` (see [Running one step](#running-one-step)).

The VQSR input for `load_vqsr` can't be downloaded: it comes from a GATK joint
call of your own samples. The step is optional, and the pipeline runs without it.

## Configure

Copy the template and fill in your own values:

```
cp gnomad_qc/v3/configuration/config.yaml my_config.yaml
```

For a filled-in starting point, use `config.example.yaml` in the same folder
instead. It is sized for a ~4,400-sample genome callset, with `step_overrides`
for each heavy step, and only the paths need changing.

| Key | What to set |
| --- | --- |
| `paths.base_dir`, `write_dir`, `log_dir`, `vds_file` | Absolute paths to your data and output locations. `code_dir` can stay empty when using the container |
| `paths.tmp_dir` | Hail's scratch directory, on node-local disk, e.g. `'${SLURM_TMPDIR}/hail'` |
| `spark.*` | Driver/executor memory and cores, matched to the resources Slurm gives the job |
| `spark.local_dir` | Spark's shuffle directory, on node-local disk, e.g. `'${SLURM_TMPDIR}/spark'` |
| `java.tool_options` | JVM options; the G1GC defaults suit large-heap Hail workloads |
| `reference.grch38_fasta`, `grch38_fasta_index` | Your copy of the GRCh38 FASTA and `.fai` |
| `reference.vep_config` | Hail VEP config JSON for the `vep` step |
| `reference.public_resource_dir` | The gnomAD resource mirror made by `host/download_resources.sh` (see [Reference tables](#reference-tables)) |
| `step_overrides` | Optional per-step Spark settings, since steps range from ~10G to ~300G peak memory |

`step_overrides` is keyed by the pipeline script's name, and any `spark` key can
be overridden:

```yaml
step_overrides:
  random_forest:        # annotate_for_rf, train_rf, apply_rf
    driver_memory: "200g"
    executor_memory: "200g"
    executor_cores: "64"
```

## Get the image

### Option 1: pull a published image (no Docker needed)

Images are published to the GitHub Container Registry and can be pulled without
a GitHub account. On the cluster:

```
module load apptainer
apptainer pull gnomad_qc_v1.0.1.sif docker://ghcr.io/canadian-genomic-data-commons/gnomad_qc_container:1.0.1
```

Use a fixed version such as `1.0.1` rather than `latest`, so every run of a
callset uses the same image. Check that it works with
`apptainer run gnomad_qc_v1.0.1.sif --help`, which lists the steps.

If the pull keeps failing partway through (`unexpected EOF`, `PROTOCOL_ERROR`),
use `host/pull_image.sh` instead. It downloads each layer with `curl`, resumes
from where it stopped when the connection drops, then builds the `.sif`:

```
module load apptainer
bash host/pull_image.sh 1.0.1      # writes gnomad_qc_v1.0.1.sif in the current directory
```

If compute nodes have no outbound network, pull on the login node and copy the
`.sif` to scratch.

Build or pull against node-local disk rather than a parallel filesystem. The
Hail and Spark layers contain tens of thousands of small files, and on Lustre a
build takes over 40 minutes instead of about 2:

```
export APPTAINER_TMPDIR=$SLURM_TMPDIR/aptmp
export APPTAINER_CACHEDIR=$SLURM_TMPDIR/apcache
```

### Option 2: build it yourself

Requires Docker, so build on a workstation or VM -- HPC login nodes do not have it.

```
docker build -t gnomad_qc:latest .
```

The build applies the fixes in `patches/gnomad_methods/` and should print
`Applied:` (or `Already applied:`) for each. On Docker Desktop, if the build
fails with `error reading from server: EOF`, give Docker's VM more memory (on
Windows, in `.wslconfig`) and retry.

Then convert for the cluster:

```
docker save gnomad_qc:latest -o gnomad_qc.tar
# transfer gnomad_qc.tar to the cluster, then, in a Slurm job:
module load apptainer
export APPTAINER_TMPDIR=$SLURM_TMPDIR/aptmp APPTAINER_CACHEDIR=$SLURM_TMPDIR/apcache
apptainer build $SLURM_TMPDIR/gnomad_qc.sif docker-archive://gnomad_qc.tar
cp $SLURM_TMPDIR/gnomad_qc.sif /path/to/scratch/
```

## Publishing a new image (maintainers)

`.github/workflows/publish-image.yml` builds the image and pushes it to ghcr.io
using the repository's own `GITHUB_TOKEN` -- no registry credentials to manage.
It does not run on every push. Trigger it either way:

```
git tag v1.0.2 && git push origin v1.0.2    # publishes 1.0.2, 1.0 and latest
```

or from the Actions tab, "Publish image" -> "Run workflow", giving a tag name.

New versions are published under the existing public package. Update the
version in the pull commands above, and the default in `host/pull_image.sh`.

## Run

Each pipeline step is selected with `--step`. Run without arguments to list the
steps in order, with the script and flags each one runs:

```
apptainer run gnomad_qc.sif
```

### Running one step

Bind your data directories into the container **at the same paths**, so the
absolute paths in your config work inside and outside the container. Include
the gnomAD resource mirror (`reference.public_resource_dir`) if it is outside
your project folder. Bind `$SLURM_TMPDIR` too if your config uses it, and pass
it through, since `--cleanenv` clears the environment:

```bash
#!/bin/bash
#SBATCH --cpus-per-task=40
#SBATCH --mem=180G
#SBATCH --time=24:00:00

module load apptainer

apptainer run --cleanenv \
    --bind /path/to/project:/path/to/project \
    --bind ${SLURM_TMPDIR}:${SLURM_TMPDIR} \
    --env GNOMAD_QC_CONFIG=/path/to/project/my_config.yaml \
    --env SLURM_TMPDIR=${SLURM_TMPDIR} \
    /path/to/gnomad_qc.sif --step sample_qc
```

Any arguments after the step name are passed to the script, e.g.
`--step compute_hard_filters --min_n_snp 10000` or `--step sample_qc --overwrite`.
Each step logs to `<script name>.log` in `log_dir`.

### Running the pipeline

Run the steps in the order `apptainer run gnomad_qc.sif` lists them; each
step reads what earlier steps wrote to `write_dir`.

### VEP

The `vep` step runs on the host, not in this container: VEP runs in its own
container, which can't be started from inside this one. It uses the same code
and config and writes `vep.ht` to `write_dir`, where `create_release_sites_ht`
reads it. Run it after `final_filter`:

```
sbatch host/run_vep.sh /path/to/my_config.yaml /path/to/python/env
```

You need:

- a Python environment on the host with Hail 0.2.134, `gnomad` (install
  `gnomad_methods` from this repository), `pyyaml` and `slackclient`
- your own VEP install and cache for GRCh38
- `reference.vep_config` pointing at a Hail VEP config JSON. Its `command`
  launches VEP (e.g. a wrapper around your VEP container), and its
  `vep_json_schema` decides which VEP fields are kept in `vep.ht`

Submit from the repository root, or set `GNOMAD_QC_REPO`. Adjust the
`#SBATCH` lines in the script to your callset.

## Resource requirements

Benchmarked on ~4400 samples on Narval/DRAC. Steps vary enormously, so size each
Slurm job individually and mirror the values in `step_overrides`.

| Step | Cores | Peak memory | Est. time (h) |
| --- | --- | --- | --- |
| *(VDS creation -- `build_vds/`, run before the container's steps)* | 40 | 300G | 120 |
| `create_last_END_positions` | 40 | 50G | 24 |
| `generate_qc_annotations` | 40 | 150G | 36 |
| `impute_sex` | 40 | 250G | 24 |
| `sample_qc` | 40 | 50G | 12 |
| `pc_relate` | 40 | 125G | 4 |
| `ancestry_classification` | 40 | 250G | 6 |
| `generate_metadata` | 40 | 10G | 0.5 |
| `find_dups` / `infer_families` / `run_mendel_errors` / `finalize_ped` | 40 | 100G | 4 |
| `freq_run_dense_annotations` / `freq_correct_high_ab_hets` | 64 | 250G | 15 |
| `annotate_for_rf` / `train_rf` / `apply_rf` | 64 | 200G | 6 |
| `create_bin_ht` / `create_aggregated_bin_ht` | 40 | 110G | 1 |
| `final_filter` | 40 | 75G | 1 |
| `vep` | 40 | 75G | 24 |
| `create_release_sites_ht` | 40 | 75G | 0.5 |

Steps not listed were not benchmarked separately. Total pipeline runtime is
roughly 7-8 days.

Other guidance from the benchmarking work:

- **Single node is preferred.** Multi-node showed poor parallelization -- a 3-node
  `compute_qc_mt` run managed only 29% CPU efficiency.
- **Use node-local storage for `tmp_dir` and `spark.local_dir`.** Lustre is slow
  at the many-small-file access patterns Hail generates.
- **Do not set `PYSPARK_SUBMIT_ARGS`** -- it overrides `hl.init()` and causes memory
  conflicts. This container drives Spark through `init_hail()` instead.
- **Resuming:** outputs of finished steps stay in `write_dir`, so after a failure
  rerun from the failed step. Pass `--overwrite` to steps that support it when
  rerunning over existing output.

## Storage

| Category | Size |
| --- | --- |
| VDS files | 5.2TB |
| Intermediate checkpoints and scratch | 10-15TB |
| QC outputs (HTs, TSVs, metrics) | ~50GB |
| Long-term archive | 4-10TB |
| **Total recommended** | **20-30TB** (minimum 10TB) |

Peak usage runs 3-4x input data size during processing.

Hail tables are stored as many small files. If your storage limits the number of
files, check that limit too: the random-forest steps alone write about 600,000
files at the default `--n_partitions` (see [Running the pipeline](#running-the-pipeline)).
