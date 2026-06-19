"""
Location of data files
======================

Use as ::

    from mda_3d_2pt.data.files import *

"""

__all__ = [
    "MDANALYSIS_LOGO",  # example file of MDAnalysis logo
]

import importlib.resources

data_directory = importlib.resources.files("mda_3d_2pt") / "data"

MDANALYSIS_LOGO = data_directory / "mda.txt"
