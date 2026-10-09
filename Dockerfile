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

# Java 11 comes from this base image. Hail 0.2.x is tested on Java 8/11; the
# Narval runs this container reproduces used Java 11 + Hail 0.2.134 + Spark 3.5.
# (The previous google/cloud-sdk:slim base rebased onto Debian trixie, which
# dropped openjdk-11 from its repos.)
FROM eclipse-temurin:11-jre-jammy

# --- System deps (Java 11 already provided by the base image) ----------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    g++ \
    less \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libncurses5-dev \
    patch \
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
ENV HAIL_VERSION="0.2.134"
# No pip cache: it only adds to the image, and to the disk the build needs.
ENV PIP_NO_CACHE_DIR=1
# Hail on its own first. It pulls in PySpark, which ships as source and has to
# be built into a wheel here -- the most memory- and disk-hungry moment of the
# build. Doing it alone keeps the peak lower, which matters on Docker Desktop,
# whose build engine crashes ("error reading from server: EOF") when its VM runs
# short. A separate layer also means a later failure doesn't redo this step.
RUN python3 -m pip install hail==${HAIL_VERSION}
RUN python3 -m pip install \
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
# Required, despite only being used for notifications: gnomad/utils/slack.py does
# `from slack import WebClient` at module level, and 62 pipeline scripts import
# from that module, so without this most steps fail at import with
# ModuleNotFoundError. Pinned to the version the pipeline was tested against.
RUN python3 -m pip install slackclient==2.5.0

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
# Fixes to gnomad_methods that aren't upstream, kept as patch files against the
# pinned submodule commit (a submodule records only a commit id, so local edits
# to it would not survive a fresh clone). A patch that is already present is
# skipped, so a build from a working tree that carries the fix still succeeds;
# one that neither applies nor is present fails the build.
COPY patches/gnomad_methods /opt/patches/gnomad_methods
RUN cd /opt/gnomad_methods && \
    for p in /opt/patches/gnomad_methods/*.patch; do \
        if patch -p1 -R --dry-run -s -f < "$p" > /dev/null 2>&1; then \
            echo "Already applied: $p"; \
        else \
            patch -p1 -s < "$p" || { echo "ERROR: could not apply $p"; exit 1; }; \
            echo "Applied: $p"; \
        fi; \
    done
RUN python3 -m pip install --no-deps -e /opt/gnomad_methods

# --- gnomad_qc pipeline scripts ---------------------------------------------
# gnomad_qc has no setup.py upstream, so it's put on PYTHONPATH instead of
# pip-installed. Note the nested destination: COPY copies the *contents* of the
# source directory, so `COPY gnomad_qc /opt/gnomad_qc` would drop v3/, v4/ and
# pipeline_config.py directly into /opt/gnomad_qc, leaving no importable
# `gnomad_qc` package and breaking the paths entrypoint.sh resolves.
COPY gnomad_qc /opt/gnomad_qc/gnomad_qc
ENV PYTHONPATH="/opt/gnomad_qc"

# --- entrypoint --------------------------------------------------------
# Dispatches `--step <name>` to the matching pipeline script. Each script
# reads paths + Spark settings from the config.yaml pointed to by
# GNOMAD_QC_CONFIG (bind-mount your own file at container run time).
COPY entrypoint.sh /opt/gnomad_qc/entrypoint.sh
RUN sed -i "s/\\r$//" /opt/gnomad_qc/entrypoint.sh && chmod +x /opt/gnomad_qc/entrypoint.sh
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
