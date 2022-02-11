import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams['axes.linewidth'] = 2
rcParams.update({'font.size': 14})
import matplotlib
matplotlib.rcParams['axes.formatter.useoffset'] = False

lys_resid = [1382, 1385, 1398, 1416,
1427, 1432, 1433, 1435,
1442, 1451, 1482, 1487,
1490
]

a = np.load('acbi1.npy') * 10
p1 = np.load('protac1.npy') * 10 
p2 = np.load('protac2.npy') * 10 

plt.figure(figsize=(20,15))
for i in range(a.shape[0]):
    plt.subplot(3,5,i+1)
    h1 = np.histogram(a[i].flatten(), bins=500, range=(0,120))
    h2 = np.histogram(p1[i].flatten(), bins=500, range=(0,120))
    h3 = np.histogram(p2[i].flatten(), bins=500, range=(0,120))

    c1 = [136/256.,169/256.,108/256.]  # smudge: acbi1
    c2 = [242/256.,186/256.,187/256.]  # salmon: protac2
    c3 = [161/256.,163/256.,198/256.]  # slate: protac1
    plt.plot(h1[1][0:-1],h1[0]/float(np.sum(h1[0])),label='ACBI1', lw=3, color=c1)
    plt.plot(h2[1][0:-1],h2[0]/float(np.sum(h2[0])),label='PROTAC 1', lw=3, color=c3)
    plt.plot(h3[1][0:-1],h3[0]/float(np.sum(h3[0])),label='PROTAC 2', lw=3, color=c2)
    plt.title('LYS%d' %lys_resid[i])
    plt.xlim(0,30)
    plt.grid(True)
    

plt.legend()
plt.grid(True)
plt.xlabel(r'Distance between LYS_N$\zeta$(SMARCA2) and C (Ub c-terminus) ($\AA$)', fontsize=12)
plt.ylabel('Probability')
plt.savefig('img_lys_density_sep_all.png',dpi=200,transparent=True)
