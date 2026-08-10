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

``import MDAnalysis as mda
import trvdos
u = mda.Universe("topol.tpr", "traj.trr")``


