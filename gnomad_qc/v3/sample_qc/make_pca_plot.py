# %%
import hail as hl
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os

from gnomad_qc.pipeline_config import BASE, WRITE_DIR

# Set plot parameters
fsz = 15
msz = 12
mrkr = 'o'
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
matplotlib.rc('font', size=fsz)

scores = os.path.join(BASE, WRITE_DIR, "ancestry_pca_scores.ht")
loadings = os.path.join(BASE, WRITE_DIR, "ancestry_pca_loadings.ht")
eigenvalues = os.path.join(BASE, WRITE_DIR, "pop_pca_eigenvalues.ht")

pca_scores = hl.read_table(scores).to_pandas()
pca_loadings = hl.read_table(loadings).to_pandas()
pca_eigenvalues = hl.read_table(eigenvalues).to_pandas()
ancestry_assignment = hl.read_table(os.path.join(BASE, WRITE_DIR, "ancestry_rf.assignment.legacy.ht")).to_pandas()

# %%
inds_BQC = [list(ancestry_assignment["s"].values).index(s) for s in list(ancestry_assignment["s"].values) if ('BQC' in s)]
inds_BCA = [list(ancestry_assignment["s"].values).index(s) for s in list(ancestry_assignment["s"].values) if (s[0] == 'F')]
inds_SIC = [list(ancestry_assignment["s"].values).index(s) for s in list(ancestry_assignment["s"].values) if not (('BQC' in s) or (s[0] == 'F'))]
# For PCA ancestry plot variant selection don't choose on the X, and have AF > 0.1%

# Extract eigenvalues and calculate percentage variance explained
eigenvalues = np.array(pca_eigenvalues['eigenvalue'])
total_variance = np.sum(eigenvalues)
percentage_variance_explained = (eigenvalues / total_variance) * 100

# Set up the PCA plot
fig_pca = plt.figure()
fig_pca.set_size_inches(10, 10)
fig_pca.set_dpi(100)

# Update axis labels with variance explained
plt.xlabel(f'PC1 ({percentage_variance_explained[0]:.2f}% variance)', fontsize=fsz)
plt.ylabel(f'PC2 ({percentage_variance_explained[1]:.2f}% variance)', fontsize=fsz)

# Legend for the three datasets
plt.scatter(-10, -10, s=msz, c=cycle[3], alpha=1, marker=mrkr, edgecolors='none', label="SickKids")
plt.scatter(-10, -10, s=msz, c=cycle[1], alpha=1, marker=mrkr, edgecolors='none', label="GSC")
plt.scatter(-10, -10, s=msz, c=cycle[0], alpha=1, marker=mrkr, edgecolors='none', label="McGill")

plt.legend(loc="lower right")

# Scatter plot for each group (PC1 vs PC2)
plt.scatter(
    [pca_scores['scores'][i][0] for i in inds_SIC],
    [pca_scores['scores'][i][1] for i in inds_SIC],
    s=msz, alpha=1.00, c=cycle[3], marker=mrkr, edgecolors='none'
)
plt.scatter(
    [pca_scores['scores'][i][0] for i in inds_BCA],
    [pca_scores['scores'][i][1] for i in inds_BCA],
    s=msz, alpha=1.00, c=cycle[1], marker=mrkr, edgecolors='none'
)
plt.scatter(
    [pca_scores['scores'][i][0] for i in inds_BQC],
    [pca_scores['scores'][i][1] for i in inds_BQC],
    s=msz, alpha=1.00, c=cycle[0], marker=mrkr, edgecolors='none'
)

# Set axis limits
ax = plt.gca()
ax.set_xlim([-0.05, 0.40])
ax.set_ylim([-0.17, 0.03])

plt.show()

# %%

# Merge ancestry assignment with PCA scores
pca_scores = pca_scores.merge(ancestry_assignment[['s', 'pop']], left_on='s', right_on='s')

# Define unique ancestry groups and assign colors
unique_ancestries = pca_scores['pop'].unique()
colors = plt.cm.get_cmap('tab10', len(unique_ancestries))

# Plot PCA with ancestry labels
fig, ax = plt.subplots(figsize=(10, 10), dpi=100)
for i, ancestry in enumerate(unique_ancestries):
    subset = pca_scores[pca_scores['pop'] == ancestry]
    ax.scatter(
        subset['scores'].apply(lambda x: x[0]), 
        subset['scores'].apply(lambda x: x[1]),
        label=ancestry,
        color=colors(i),
        alpha=0.7,
        edgecolors='none'
    )

# Add labels with explained variance
plt.xlabel(f'PC1 ({percentage_variance_explained[0]:.2f}% variance)', fontsize=fsz)
plt.ylabel(f'PC2 ({percentage_variance_explained[1]:.2f}% variance)', fontsize=fsz)

# Add legend
ax.legend(title="Ancestry")

# Set axis limits
ax = plt.gca()
ax.set_xlim([-0.05, 0.40])
ax.set_ylim([-0.17, 0.03])

plt.show()
fig_pca.savefig(os.path.join(BASE, WRITE_DIR, "PCA_by_dataset.png"))
fig.savefig(os.path.join(BASE, WRITE_DIR, "PCA_by_ancestry.png"))

# Zoom in
ax.set_xlim([-0.05, 0.075])
fig.savefig(os.path.join(BASE, WRITE_DIR, "PCA_by_ancestry_zoomed.png"))