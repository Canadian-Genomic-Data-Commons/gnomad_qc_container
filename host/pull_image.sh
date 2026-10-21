#!/bin/bash
# Downloads the pipeline image from ghcr.io and builds a .sif from it. Use this
# instead of `apptainer pull` when the pull keeps failing partway through
# ("unexpected EOF", "PROTOCOL_ERROR"): each layer is fetched with curl and
# resumed from where it stopped if the connection drops.
#
# Usage:
#   bash host/pull_image.sh [version] [output.sif]
#   (defaults: 1.0.1 and gnomad_qc_v<version>.sif in the current directory)
set -euo pipefail
shopt -s inherit_errexit

TAG="${1:-1.0.1}"
OUT="${2:-gnomad_qc_v${TAG}.sif}"
REPO="${GHCR_REPO:-canadian-genomic-data-commons/gnomad_qc_container}"
REGISTRY="https://ghcr.io"
LAYOUT="${OUT%.sif}.oci"          # downloaded files are kept here until the build succeeds
MAX_TRIES=30

mkdir -p "$LAYOUT/blobs/sha256"

get_token() {
    local url="$REGISTRY/token?scope=repository:$REPO:pull&service=ghcr.io"
    local resp
    if [[ -n "${APPTAINER_DOCKER_PASSWORD:-}" ]]; then
        resp=$(curl -fsS -u "${APPTAINER_DOCKER_USERNAME:?set APPTAINER_DOCKER_USERNAME}:$APPTAINER_DOCKER_PASSWORD" "$url")
    else
        resp=$(curl -fsS "$url")
    fi || { echo "ERROR: ghcr.io refused the request. Check the version exists, and for a private copy of the image set APPTAINER_DOCKER_USERNAME and APPTAINER_DOCKER_PASSWORD." >&2; exit 1; }
    python3 -I -c 'import json,sys; print(json.load(sys.stdin)["token"])' <<<"$resp"
}

ACCEPT="application/vnd.oci.image.index.v1+json,application/vnd.oci.image.manifest.v1+json,application/vnd.docker.distribution.manifest.list.v2+json,application/vnd.docker.distribution.manifest.v2+json"

# Writes a manifest (or index) to the layout and prints its digest.
fetch_manifest() {
    local ref="$1" tmp="$LAYOUT/manifest.tmp" digest
    curl -fsS -H "Authorization: Bearer $TOKEN" -H "Accept: $ACCEPT" \
        "$REGISTRY/v2/$REPO/manifests/$ref" -o "$tmp"
    digest="sha256:$(sha256sum "$tmp" | cut -d' ' -f1)"
    mv "$tmp" "$LAYOUT/blobs/sha256/${digest#sha256:}"
    echo "$digest"
}

# Downloads one blob, resuming after a dropped connection, and checks its digest.
fetch_blob() {
    local digest="$1" size="$2" file="$LAYOUT/blobs/sha256/${1#sha256:}" try have
    for ((try = 1; try <= MAX_TRIES; try++)); do
        have=$(stat -c %s "$file" 2>/dev/null || echo 0)
        if (( have == size )); then
            if [[ "$(sha256sum "$file" | cut -d' ' -f1)" == "${digest#sha256:}" ]]; then
                return 0
            fi
            echo "  digest mismatch for $digest, downloading again" >&2
            rm -f "$file"; have=0
        elif (( have > size )); then
            rm -f "$file"; have=0
        fi
        (( try > 1 )) && echo "  resuming $digest at $((have / 1048576)) of $((size / 1048576)) MiB (try $try)" >&2
        TOKEN=$(get_token)
        curl -fsSL -C - --retry 3 --retry-delay 5 \
            -H "Authorization: Bearer $TOKEN" \
            "$REGISTRY/v2/$REPO/blobs/$digest" -o "$file" || sleep 5
    done
    echo "ERROR: could not download $digest after $MAX_TRIES tries" >&2
    exit 1
}

TOKEN=$(get_token)
top=$(fetch_manifest "$TAG")

# CI publishes an index (amd64 image plus build attestations); pick the amd64 image.
read -r media manifest <<<"$(python3 -I - "$LAYOUT/blobs/sha256/${top#sha256:}" "$top" <<'EOF'
import json, sys
path, digest = sys.argv[1], sys.argv[2]
doc = json.load(open(path))
media = doc.get("mediaType", "application/vnd.oci.image.manifest.v1+json")
if "manifests" in doc:
    for m in doc["manifests"]:
        p = m.get("platform", {})
        if p.get("os") == "linux" and p.get("architecture") == "amd64":
            print(m["mediaType"], m["digest"]); break
    else:
        sys.exit("no linux/amd64 image in the index")
else:
    print(media, digest)
EOF
)"
[[ "$manifest" != "$top" ]] && fetch_manifest "$manifest" >/dev/null

# Config and layers listed in the image manifest, as "digest size" lines.
mapfile -t blobs < <(python3 -I -c '
import json, sys
m = json.load(open(sys.argv[1]))
for b in [m["config"]] + m["layers"]:
    print(b["digest"], b["size"])
' "$LAYOUT/blobs/sha256/${manifest#sha256:}")

n=0
for line in "${blobs[@]}"; do
    n=$((n + 1))
    read -r digest size <<<"$line"
    echo "[$n/${#blobs[@]}] $digest ($((size / 1048576)) MiB)"
    fetch_blob "$digest" "$size"
done

manifest_size=$(stat -c %s "$LAYOUT/blobs/sha256/${manifest#sha256:}")
echo '{"imageLayoutVersion":"1.0.0"}' > "$LAYOUT/oci-layout"
printf '{"schemaVersion":2,"manifests":[{"mediaType":"%s","digest":"%s","size":%s}]}\n' \
    "$media" "$manifest" "$manifest_size" > "$LAYOUT/index.json"

echo "All layers downloaded. Building $OUT ..."
apptainer build --force "$OUT" "oci:$LAYOUT"
rm -rf "$LAYOUT"
echo "Done: $OUT"
