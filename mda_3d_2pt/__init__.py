"""
MDA-3D-2PT
Development of 3D-2PT package for MDAnalysis (with compatibility for complex solvent molecules)
"""

# Add imports here
from importlib.metadata import version

__version__ = version("mda_3d_2pt")
from .vdos import *
