# Dockerfile for running the gnomad_qc pipeline scripts (v3/v4/v5) plus the
# gnomad_methods ("gnomad") library outside of Dataproc, e.g. on an HPC
# cluster via Apptainer/Singularity.
#
# Build context must be the repo root that contains both `gnomad_methods/`
# and `gnomad_qc/` as sibling directories (this file's directory).
#
#   docker build -t gnomad_qc:latest .
#
# See the bottom of this file for the Apptainer/Singularity conversion steps.

FROM google/cloud-sdk:slim

# --- System deps + Java 11 (required by Hail/Spark) -------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    g++ \
    less \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libncurses5-dev \
    openjdk-11-jre-headless \
    pkg-config \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    unzip \
    wget \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --upgrade pip

# --- Python deps --------------------------------------------------------
# Pin the Hail version to whatever you use to read/write your existing
# Hail Tables/MatrixTables/VDS - Hail's on-disk format is version-sensitive.
ENV HAIL_VERSION="0.2.113"
RUN python3 -m pip install \
    hail==${HAIL_VERSION} \
    numpy \
    pandas \
    scipy \
    scikit-learn \
    matplotlib \
    seaborn \
    networkx \
    onnx \
    onnxruntime \
    skl2onnx \
    annoy \
    bokeh \
    tabulate \
    pyyaml \
    jsonschema \
    papermill \
    ipython \
    click
# Optional: Slack notifications used by gnomad.utils.slack (`from slack import
# WebClient`) - only needed if you actually use that helper.
# RUN python3 -m pip install slackclient

# --- GCS connector (only needed if you read/write gs:// paths from here) -----
RUN export SPARK_HOME=$(find_spark_home.py) && \
    curl https://storage.googleapis.com/hadoop-lib/gcs/gcs-connector-hadoop2-2.0.1.jar \
         >$SPARK_HOME/jars/gcs-connector-hadoop2-2.0.1.jar && \
    mkdir -p $SPARK_HOME/conf && \
    touch $SPARK_HOME/conf/spark-defaults.conf && \
    sed -i $SPARK_HOME/conf/spark-defaults.conf \
        -e 's:spark\.hadoop\.google\.cloud\.auth\.service\.account\.enable.*:spark.hadoop.google.cloud.auth.service.account.enable true:' \
        -e 's:spark\.hadoop\.google\.cloud\.auth\.service\.account\.json\.keyfile.*:spark\.hadoop\.google\.cloud\.auth\.service\.account\.json\.keyfile /gsa-key/key.json:'

# --- gnomad_methods (the "gnomad" package that gnomad_qc imports) -----------
# This is a git submodule. A clone without `--recursive` leaves it empty, which
# would otherwise build a broken image that only fails at pipeline run time.
COPY gnomad_methods /opt/gnomad_methods
RUN test -f /opt/gnomad_methods/setup.py || ( \
        echo "ERROR: gnomad_methods/ is empty. Run: git submodule update --init" && \
        exit 1 )
RUN python3 -m pip install --no-deps -e /opt/gnomad_methods

# --- gnomad_qc pipeline scripts ---------------------------------------------
# gnomad_qc has no setup.py upstream, so it's put on PYTHONPATH instead of
# pip-installed.
COPY gnomad_qc /opt/gnomad_qc
ENV PYTHONPATH="/opt/gnomad_qc:${PYTHONPATH}"

# --- entrypoint --------------------------------------------------------
# Dispatches `--step <name>` to the matching pipeline script. Each script
# reads paths + Spark settings from the config.yaml pointed to by
# GNOMAD_QC_CONFIG (bind-mount your own file at container run time).
COPY entrypoint.sh /opt/gnomad_qc/entrypoint.sh
RUN chmod +x /opt/gnomad_qc/entrypoint.sh
ENTRYPOINT ["/opt/gnomad_qc/entrypoint.sh"]

WORKDIR /home

# ---------------------------------------------------------------------------
# Convert to Apptainer/Singularity (run these on a machine with Docker, e.g.
# your laptop or a cloud VM - most HPC login/compute nodes don't have Docker):
#
#   docker build -t gnomad_qc:latest .
#   docker push <your-registry>/gnomad_qc:latest   # e.g. Docker Hub, ghcr.io
#
# Then on the cluster (DRAC/Compute Canada nodes ship "apptainer", which is
# singularity-compatible - use the `apptainer` command, or `module load
# apptainer` first):
#
#   apptainer pull gnomad_qc.sif docker://<your-registry>/gnomad_qc:latest
#
# If the compute nodes have no internet access but the login node does, pull
# there and copy the .sif to scratch. If you don't want to use a registry at
# all, save+transfer the image instead:
#
#   docker save gnomad_qc:latest -o gnomad_qc.tar
#   # scp gnomad_qc.tar to the cluster
#   apptainer build gnomad_qc.sif docker-archive://gnomad_qc.tar
#
# Running a step from the .sif: bind your data dir and your filled-in
# config.yaml, and point GNOMAD_QC_CONFIG at where you bound the config to.
#
#   apptainer run \
#       --bind /lustre07/scratch/kobe/gnomad:/data \
#       --bind /path/to/my_config.yaml:/config/config.yaml \
#       --env GNOMAD_QC_CONFIG=/config/config.yaml \
#       gnomad_qc.sif --step ancestry_classification
#
# See gnomad_qc/v3/configuration/config.yaml for the config template (paths +
# Spark settings) and entrypoint.sh for the full list of available steps.
# ---------------------------------------------------------------------------
