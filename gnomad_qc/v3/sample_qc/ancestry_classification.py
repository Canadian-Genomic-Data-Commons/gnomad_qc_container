# %% 
"""
Use random forest (RF) model on PCs for samples with known ancestry labels from gnomad dataset and apply that RF model to assign ancestry to HostSeq samples. 
Could alternatively train our own model on HostSeq samples with known ancestry labels assuming we have enough of those in the clinical records.
"""
import os

import onnx
import hail as hl
from gnomad.sample_qc.ancestry import apply_onnx_classification_model, assign_population_pcs
from gnomad.utils.filtering import filter_to_adj

from gnomad_qc.pipeline_config import BASE, WRITE_DIR, VDS_FILE, init_hail

read_if_exists = True

# Configure Hail using the shared config.yaml (paths + Spark settings).
init_hail("ancestry_classification")
# %%
# v3.1 PCA loadings.
gnomad_v3_loadings = os.path.join(BASE, "hail", "gnomad.v3.1.pca_loadings.ht")
# v3.1 ONNX RF model.
gnomad_v3_onnx_rf = os.path.join(BASE, "hail", "gnomad.v3.1.RF_fit.onnx")
# Test dataset to apply projection and genetic ancestry group assignment to
vds_to_project = os.path.join(BASE, VDS_FILE)

# v3.1 output paths.
test_mt_output_path = os.path.join(BASE, WRITE_DIR, "gnomad_v3.1_ancestry_rf.mt")
test_scores_output_path = os.path.join(BASE, WRITE_DIR, "gnomad_v3.1_ancestry_rf.scores.ht")
gnomad_v3_assignment_path = os.path.join(BASE, WRITE_DIR, "gnomad_v3.1_ancestry_rf.assignment.ht")

# Define the number of PCs used for v2 and v3 genetic ancestry group classification
v3_num_pcs = 16
# Define the RF minimum probability used for v2 and v3 genetic ancestry group classification
v3_min_prob = 0.75

# %% Load ONNX models
print("Loading gnomAD v3.1 ONNX RF model")
with hl.hadoop_open(gnomad_v3_onnx_rf, "rb") as f:
    v3_onx_fit = onnx.load(f)

# %% Load gnomAD v3.1 loadings Hail Table and the VariantDataset to apply projection and genetic ancestry group assignment to
print("Loading gnomAD v3.1 loadings Hail Table and VariantDataset")
vds = hl.vds.read_vds(vds_to_project)
v3_loading_ht = hl.read_table(gnomad_v3_loadings)

# %% Perform PC projection using the v3.1 PCA loadings 

# Create dense MatrixTable of only the variants found in the loadings Table
# Reduce variant data to only annotations needed for split and densify.
# This includes annotations needed for our standard genotype filter ('adj').
vds = hl.vds.VariantDataset(
    vds.reference_data, 
    vds.variant_data.select_entries("LA", "LGT", "GQ", "DP", "LAD")
)

# Split multiallelics.
vds = hl.vds.split_multi(vds, filter_changed_loci=True)

# Filter to variants in the loadings Table.
vds = hl.vds.filter_variants(vds, v3_loading_ht)

# Densify VDS.
mt = hl.vds.to_dense_mt(vds)

# %% We recommend filtering to entries meeting GQ, DP and het AB 'adj' thresholds
mt = filter_to_adj(mt)

# %% Checkpoint dense MT for PC projection
print("Checkpointing dense MT for PC projection")
mt = mt.checkpoint(
    test_mt_output_path, 
    overwrite=not read_if_exists, 
    _read_if_exists=read_if_exists
)

# %% Project test dataset genotypes onto gnomAD v3.1 loadings and checkpoint the scores
# Project new genotypes onto loadings.
print("Projecting test dataset genotypes onto gnomAD v3.1 loadings")
v3_pcs_ht = hl.experimental.pc_project(
    mt.GT, v3_loading_ht.loadings, v3_loading_ht.pca_af,
)

# Checkpoint PC projection results.
print("Checkpointing PC projection results")
v3_pcs_ht = v3_pcs_ht.checkpoint(
    test_scores_output_path, 
    overwrite=not read_if_exists, 
    _read_if_exists=read_if_exists
)

# %% Assign genetic ancestry group using ONNX RF model v3.1 RF model
print("Assigning genetic ancestry group using ONNX RF model")
ht, model = assign_population_pcs(
    v3_pcs_ht,
    pc_cols=v3_pcs_ht.scores[:v3_num_pcs],
    fit=v3_onx_fit,
    min_prob=v3_min_prob,
    apply_model_func = apply_onnx_classification_model,
)
print("Writing genetic ancestry group assignments to Hail Table")
ht = ht.checkpoint(
    gnomad_v3_assignment_path, 
    overwrite=not read_if_exists, 
    _read_if_exists=read_if_exists
)
ht.aggregate(hl.agg.counter(ht.pop))

print(f"Reading HT from {gnomad_v3_assignment_path}")
ht = hl.read_table(gnomad_v3_assignment_path)
ht = ht.select_globals()
bgz_write = os.path.join(
    BASE, 
    WRITE_DIR, 
    "ancestry_rf.assignment.tsv.bgz"
)
print(f"Writing TSV to {bgz_write}")
ht.export(bgz_write)