Getting Started
===============

This page details how to get started with trvdos. 

Installation
************

trvdos is available on PyPI, and can be installed via

``pip install trvdos``

trvdos requires the GNU Scientific Library (GSL). Since pip cannot automatically handle C prerequisites, GSL also needs to be installed to successfully run trvdos.

GSL can be installed and managed via conda / mamba:

``mamba install gsl``


Example Usage
*************

trvdos uses MDAnalysis AtomGroups to handle MD trajectory data. As an example, we'll set up an MDAnalysis Universe using MD data stored in ``topol.tpr`` and ``traj.trr``:

.. code-block:: python

   import MDAnalysis as mda
   import trvdos
   u = mda.Universe("topol.tpr", "traj.trr")

trvdos additionally requires that MD trajectory data is unwrapped (made whole) across periodic boundary conditions. trvdos includes a helper function trvdos.unwrap to do this if necessary.

``unwrap = trvdos.unwrap(u)

Here we'll select the POPC lipid molecules for VDoS analysis and create a vdos object with an autocorrelation lag time (in number of frames) of 200:

.. code-block:: python

   sel = u.select_atoms('resname POPC')
   vdos = trvdos.vdos(sel, 200)


We can now iterate over the trajectory frame-by-frame to calculate single-frame values:

.. code-block:: python

   for ts in (u.trajectory):
     unwrap.single_frame()
     vdos.single_frame(tStep,ts.time)
     tStep += 1

The residue properties, VaCF and VDoS can then be accessed in the VDoS object and/or output to data files:

.. code-block:: python

   vdos.copyResidueList()
   vdos.postProcess(vdos.residueListCopy)
   vdos.outputGeometry("residueProperties.dat",vdos.residueListCopy)
   vdos.outputVACF("VACF.dat")
   vdos.outputVDoS("VDoS.dat")

