"""
trvdos
Development of 3D-2PT package for MDAnalysis (with compatibility for complex solvent molecules)
"""

# Add imports here
from importlib.metadata import version

__version__ = version("trvdos")
from .vdos import *
