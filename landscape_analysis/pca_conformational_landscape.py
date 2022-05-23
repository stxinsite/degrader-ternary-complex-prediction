import mdtraj as md
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import pickle, csv


# compute distances
print("Loading trajectory")
traj = md.load('cluster.mol-center-lig.xtc', top='0ps.pdb')
nres1 = len(traj.top.to_fasta()[0])
nres2 = len(traj.top.to_fasta()[1])
res1 = np.arange(nres1)
res2 = np.arange(nres1, nres1+nres2)
pairs = [[r1, r2] for r1 in res1 for r2 in res2]
print("Computing distances")
dist, pairs = md.compute_contacts(traj, pairs, scheme='closest-heavy')
is_contact = np.sum(dist < 0.5, axis=0)>0
interface_pairs = pairs[is_contact]
interface_dist = dist[:, is_contact]

# perform pca decomposition on interface distances
print("Performing pca decomposition")
pca=PCA(n_components=20)
feat_traj=pca.fit_transform(interface_dist)
nfeat = np.where(np.cumsum(pca.explained_variance_ratio_) > 0.95)[0][0]+1

# kmeans with a pre-determined number of clusters
print("Performing KMeans clustering")
km=KMeans(n_clusters=6)
km.fit(feat_traj[:, :nfeat])
km.cluster_centers_

# plot landscape
fig, ax= plt.subplots(nfeat-1, nfeat-1, figsize=(10, 10))
for iax in range(nfeat-1):
    for jax in range(iax+1, nfeat):
        ax[iax, jax-1].hist2d(*feat_traj[:, [iax, jax]].T, bins=100, 
                density=True, norm=mpl.colors.LogNorm())
        ax[iax, jax-1].scatter(*km.cluster_centers_[:, [iax, jax]].T, 
                marker='x', s=20, c='Red')
plt.tight_layout()
fig.savefig('clusters.landscape.png')
plt.close()

# identify frames closest to cluster centers
idx_center = np.argmin(
        np.sum((np.expand_dims(feat_traj, 1) - 
                np.expand_dims(km.cluster_centers_, 0))**2, axis=-1), 
        axis=0)
for iframe, frame in enumerate(idx_center):
    trj[frame].save_pdb(f"cluster_center.{iframe}.pdb")
