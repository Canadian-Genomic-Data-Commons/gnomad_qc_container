# gnomad_qc_container

Containerized [gnomAD QC pipeline](https://github.com/broadinstitute/gnomad_qc) for
running on HPC clusters (Slurm + Apptainer/Singularity) instead of Google Dataproc.

Derived from the Broad Institute's `gnomad_qc`, with the `gnomad_methods` library
included as a pinned submodule. Site-specific paths and Spark settings have been
factored out of the pipeline scripts into a single `config.yaml`, so each site
points the container at its own storage and cluster sizing without editing code.

The pipeline starts from a Hail VDS of your callset (GRCh38) and ends with a
release sites table, `release_sites.ht`.

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

## Inputs

| Input | Where it goes |
| --- | --- |
| Hail VDS of your callset (GRCh38) | `paths.vds_file` in the config |
| GRCh38 FASTA and `.fai` | `reference.grch38_fasta` / `grch38_fasta_index` |
| Sample metadata TSV, named `meta_test.tsv` | Your `write_dir`, before the first step |
| VEP install and cache | Used by the `vep` step on the host (see [VEP](#vep)) |

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

## Configure

Copy the template and fill in your own values:

```
cp gnomad_qc/v3/configuration/config.yaml my_config.yaml
```

| Key | What to set |
| --- | --- |
| `paths.base_dir`, `write_dir`, `log_dir`, `vds_file` | Absolute paths to your data and output locations. `code_dir` can stay empty when using the container |
| `paths.tmp_dir` | Hail's scratch directory, on node-local disk, e.g. `'${SLURM_TMPDIR}/hail'` |
| `spark.*` | Driver/executor memory and cores, matched to the resources Slurm gives the job |
| `spark.local_dir` | Spark's shuffle directory, on node-local disk, e.g. `'${SLURM_TMPDIR}/spark'` |
| `java.tool_options` | JVM options; the G1GC defaults suit large-heap Hail workloads |
| `reference.grch38_fasta`, `grch38_fasta_index` | Your copy of the GRCh38 FASTA and `.fai` |
| `reference.vep_config` | Hail VEP config JSON for the `vep` step |
| `reference.public_resource_dir` | Optional local copy of gnomAD's public resources, for compute nodes without internet |
| `step_overrides` | Optional per-step Spark settings, since steps range from ~10G to ~300G peak memory |

Environment variables in the config are expanded when a step starts, so
`${SLURM_TMPDIR}` resolves to each job's own node-local directory. Both
`tmp_dir` and `local_dir` should be node-local: on a parallel filesystem such as
Lustre, shuffle-heavy steps are very slow.

`step_overrides` is keyed by the pipeline script's name, and any `spark` key can
be overridden:

```yaml
step_overrides:
  random_forest:        # annotate_for_rf, train_rf, apply_rf
    driver_memory: "200g"
    executor_memory: "200g"
    executor_cores: "64"
```

The script names are `sample_qc`, `generate_qc_annotations`,
`create_last_END_positions`, `ancestry_classification`, `create_fam`,
`generate_freq`, `random_forest`, `evaluation`, `final_filter` and
`create_release_sites_ht`. `apptainer run gnomad_qc.sif` shows which script
each step runs.

## Get the image

### Option 1: pull a published image (no Docker needed)

Images are published to the GitHub Container Registry. On the cluster:

```
module load apptainer
apptainer pull gnomad_qc.sif docker://ghcr.io/canadian-genomic-data-commons/gnomad_qc_container:latest
```

If the package is private, authenticate first with a GitHub personal access token
that has `read:packages`:

```
apptainer registry login --username <github-username> docker://ghcr.io
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
git tag v1.0.0 && git push origin v1.0.0    # publishes 1.0.0, 1.0 and latest
```

or from the Actions tab, "Publish image" -> "Run workflow", giving a tag name.

The package inherits the repository's visibility, so on a private repository it
starts private. To let other labs pull it without a token, change it under
Packages -> gnomad_qc_container -> Package settings -> Change visibility.

## Run

Each pipeline step is selected with `--step`. Run without arguments to list the
steps in order, with the script and flags each one runs:

```
apptainer run gnomad_qc.sif
```

### Running one step

Bind your data directories into the container **at the same paths**, so the
absolute paths in your config work inside and outside the container. Bind
`$SLURM_TMPDIR` too if your config uses it, and pass it through, since
`--cleanenv` clears the environment:

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
step reads what earlier steps wrote to `write_dir`. Steps that need extra
arguments from you:

| Step | Arguments |
| --- | --- |
| `compute_hard_filters` | `--min_n_snp <n>` for an exome, panel or subset. The defaults are genome-wide counts (minimum 500,000 SNVs per sample), and below them every sample is filtered out |
| `apply_rf`, `create_bin_ht`, `create_aggregated_bin_ht` | `--model_id <id>`. `train_rf` writes the id to `write_dir/rf_runs/latest_model_id` |
| `final_filter` | `--model_id <id> --model_name RF --score_name RF --snp_bin_cutoff <n> --indel_bin_cutoff <n>`, e.g. cutoffs 90 and 80 |
| `load_vqsr` | `--vqsr_vcf_path <path> --header_path <path>`. Optional; skip it if you have no VQSR VCF |

`generate_metadata` flags sample outliers with the regression model from
`apply_regressed_filters` by default. To use the population-stratified model
from `apply_stratified_filters` instead, run `generate_metadata`, then
`apply_stratified_filters`, then `generate_metadata` again with
`--stratified_metrics_outlier` (the stratified filters need the populations the
first `generate_metadata` writes).

Several flags of the same script can be combined in one call, e.g.
`--step sample_qc --impute_sex --compute_hard_filters`. They run in the
script's own fixed order, so only combine steps that are adjacent in the list
above.

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
