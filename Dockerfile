# Image for running the gnomad_qc pipeline and the gnomad_methods ("gnomad")
# library on an HPC cluster via Apptainer/Singularity.
#
# Build from the repo root (this file's directory), after cloning with
# --recursive so gnomad_methods/ is populated.
#
#   docker build -t gnomad_qc:latest .
#
# See the bottom of this file for the Apptainer/Singularity conversion steps.

# Java 11 (required by Hail) comes from the base image.
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
ENV PIP_NO_CACHE_DIR=1
# Hail is installed on its own because building PySpark is the most memory-
# hungry part of the build. If Docker Desktop fails here with "error reading
# from server: EOF", give its VM more memory (on Windows, in .wslconfig).
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
# Required even though Slack notifications aren't used: most pipeline scripts
# import gnomad's Slack helper, which fails to import without this package.
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
# Fails the build if the submodule is empty (repo cloned without --recursive).
COPY gnomad_methods /opt/gnomad_methods
RUN test -f /opt/gnomad_methods/setup.py || ( \
        echo "ERROR: gnomad_methods/ is empty. Run: git submodule update --init" && \
        exit 1 )
# Applies the fixes in patches/gnomad_methods/ to the pinned gnomad_methods
# commit. Patches already applied are skipped; a patch that can't be applied
# fails the build. Windows line endings are stripped first so a build from a
# Windows checkout behaves the same as one from Linux.
COPY patches/gnomad_methods /opt/patches/gnomad_methods
RUN find /opt/gnomad_methods -type f -name '*.py' -exec sed -i 's/\r$//' {} + && \
    cd /opt/gnomad_methods && \
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
# gnomad_qc isn't pip-installable, so it's put on PYTHONPATH instead.
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
#       --bind /path/to/your/data:/data \
#       --bind /path/to/my_config.yaml:/config/config.yaml \
#       --env GNOMAD_QC_CONFIG=/config/config.yaml \
#       gnomad_qc.sif --step ancestry_classification
#
# See gnomad_qc/v3/configuration/config.yaml for the config template (paths +
# Spark settings) and entrypoint.sh for the full list of available steps.
# ---------------------------------------------------------------------------
