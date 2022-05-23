import mdtraj as md
import numpy as np
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt

def cal(lys_trr,lys_ref_pdb,ub_pdb,outputname):
    print(name)
    t = md.load(lys_trr,top=lys_ref_pdb)
    u = md.load(ub_pdb)

    nz = t.topology.select('name NZ')
    nn = t.xyz[:,nz,:]
    print(nn.shape)
    nn = nn[::10]
    uu = u.xyz[:,0,:]
    print(nn.shape, uu.shape)
    dd = []
    for i in range(1,14):
        print(i)
        try: dd.append(distance_matrix(nn[:,i,:],uu))
        except: print("skipped %d" %i)
    d3 = np.array(dd)
    np.save(outputname,d3)


cal('acbi1.trr','acbi1.pdb','ub.pdb','lys_dist.npy')
