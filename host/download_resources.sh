#!/bin/bash
# Downloads the gnomAD reference tables the pipeline reads, from gnomAD's free
# public bucket (gs://gcp-public-data--gnomad; no Google account needed), and
# lays them out where the pipeline looks for them:
#
#   <base_dir>/hail/   tables the pipeline scripts read
#   <base_dir>/        gnomad.genomes.r3.0.sites.ht and purcell5k.ht
#   <mirror_dir>/      links in gnomAD's bucket layout, for the gnomad library
#                      (set reference.public_resource_dir to this folder)
#
# Needs gsutil (part of the Google Cloud SDK). Run it on a node with internet
# access, e.g. a login node, with nohup since the downloads take hours:
#
#   nohup bash host/download_resources.sh /path/to/base_dir /path/to/gnomad_public_mirror > download_resources.log 2>&1 &
#
# Space: about 360 GB, almost all of it gnomad.genomes.r3.0.sites.ht (~340 GB)
# and dbSNP (~19 GB). If it is interrupted, rerun it: only missing files are
# downloaded.
set -euo pipefail

BASE_DIR=${1:?usage: download_resources.sh /path/to/base_dir /path/to/gnomad_public_mirror}
MIRROR_DIR=${2:?usage: download_resources.sh /path/to/base_dir /path/to/gnomad_public_mirror}
BUCKET=gs://gcp-public-data--gnomad

# Copies one table (a folder) or file from the bucket into dest_dir. rsync only
# copies what is missing, so a rerun finishes an interrupted download.
download() {
    local src="$1" dest_dir="$2" name
    name=$(basename "$src")
    echo "Downloading $src"
    mkdir -p "$dest_dir"
    if [[ "$name" == *.ht ]]; then
        gsutil -o "GSUtil:parallel_thread_count=1" -m -q rsync -r "$BUCKET/$src" "$dest_dir/$name"
    else
        gsutil -q cp -n "$BUCKET/$src" "$dest_dir/"
    fi
}

# Read by the pipeline scripts from base_dir/hail.
HAIL_TABLES=(
    resources/grch38/telomeres_and_centromeres/hg38.telomeresAndMergedCentromeres.ht
    resources/grch38/lcr_intervals/LCRFromHengHg38.ht
    resources/grch38/seg_dup_intervals/GRCh38_segdups.ht
    resources/grch38/clinvar/clinvar_20190923.ht
    resources/grch38/dbsnp/dbsnp_b151_grch38_all_20180418.ht
    resources/grch38/truth_sites.ht
    release/3.1/pca/gnomad.v3.1.pca_loadings.ht
    release/3.1/pca/gnomad.v3.1.RF_fit.onnx
)
for t in "${HAIL_TABLES[@]}"; do
    download "$t" "$BASE_DIR/hail"
done

# Read by the pipeline scripts from base_dir.
download resources/grch38/purcell_5k_intervals/purcell5k.ht "$BASE_DIR"
download release/3.0/ht/genomes/gnomad.genomes.r3.0.sites.ht "$BASE_DIR"

# The gnomad library looks for these at the same paths as in gnomAD's bucket,
# under reference.public_resource_dir. Link them rather than download twice.
link() {
    local bucket_path="$1" target="$2"
    mkdir -p "$MIRROR_DIR/$(dirname "$bucket_path")"
    ln -sfn "$target" "$MIRROR_DIR/$bucket_path"
}
for t in "${HAIL_TABLES[@]:0:5}"; do
    link "$t" "$BASE_DIR/hail/$(basename "$t")"
done
link resources/grch38/purcell_5k_intervals/purcell5k.ht "$BASE_DIR/purcell5k.ht"

echo "Done. Set reference.public_resource_dir to $MIRROR_DIR, and bind both"
echo "$BASE_DIR and $MIRROR_DIR into the container at the same paths."
