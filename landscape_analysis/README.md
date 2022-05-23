# PCA Analysis of conformational landscape

1. Perform HREMD simulations of desired system
2. Place trajectory and topology file containing only solute in this directory, with the trajectory named `clusters.mol-center-lig.xtc` and the topology names `0ps.top`. 
3. Run the python script `pca_conformational_landscape.py`. The number of clusters used in the KMeans algorithm may need to be updated so that the cluster centers coincide nicely with the local minima of the free energy landscape. 

Running the script generates PDB files of the structures of the cluster centers and a figure showing the conformational landscape projected onto all pairs of axes.

This method was used to identify local minima of the free energy landscape from HREMD simulations for figures 7, S14, and S15 
