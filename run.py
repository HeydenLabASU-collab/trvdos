# %%
import ctypes as ct
import MDAnalysis as mda
import MDA_unwrap_PBC as pbc
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import vdos as vd

# %%
TOPOL = "/Users/matthiasheyden/Dropbox (ASU)/ASU-Research/POPC/run-NVE.tpr"
TRAJ = "/Users/matthiasheyden/Dropbox (ASU)/ASU-Research/POPC/run-NVE_test.trr"
u = mda.Universe(TOPOL,TRAJ)
pbc = pbc.unwrap(u)
sel = u.select_atoms("resname POPC")
vdos = vd.vdos(sel,200)

# %%
vd.vdosLib.omp_set_num_threads(4)
tStep = 0
for ts in tqdm(u.trajectory):
    pbc.run()
    vdos.processStep(tStep,ts.time)
    tStep += 1

# %%
vdos.postProcess()
vdos.outputGeometry("residueProperties.dat")
vdos.outputVACF("VACF.dat")
vdos.outputVDoS("VDoS.dat")

# %%
# plot
data=np.array(vdos.totVDoS[0])
fig, ax = plt.subplots()
ax.plot(vdos.tau, data, linewidth = 2.0)
ax.set(xlim=(0, 1.8), xticks=np.arange(0,1.8,0.2),
       ylim=(0, 400000.0), yticks=np.arange(0,400001,100000))
plt.show()


