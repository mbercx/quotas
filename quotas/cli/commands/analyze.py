# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os
import matplotlib.pyplot as plt

from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Locpot, Vasprun, Outcar
from pymatgen.analysis.surface_analysis import WorkFunctionAnalyzer

"""
Module that defines the analysis tools for the quotas package.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.2"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Apr 2018"


def wf(directory, plot_potential=False):

    directory = os.path.abspath(directory)

    wf_analyzer = WorkFunctionAnalyzer.from_files(
        poscar_filename=os.path.join(directory, "POSCAR"),
        locpot_filename=os.path.join(directory, "LOCPOT"),
        outcar_filename=os.path.join(directory, "OUTCAR")
    )

    if plot_potential:
        wf_analyzer.get_locpot_along_slab_plot().show()
        pass #TODO Finish this part

    work_function = wf_analyzer.work_function
    print("Work function = " + str(work_function) + " eV")


def bandgap(directory):

    # Load the vasprun.xml file
    vasprun = Vasprun(os.path.join(directory, "vasprun.xml"))

    # Extract the band gap from the band structure
    band_gap = vasprun.get_band_structure().get_band_gap()

    print("Band gap = " + str(band_gap) + " eV")

def pdos(vasprun):

    # Load the vasprun.xml file
    out = Vasprun(os.path.abspath(vasprun))

    # Extract the types of species in the structure
    species = out.initial_structure.types_of_specie

    # Extract the complete DOS
    element_dos = out.complete_dos.get_element_dos()

    # Initialize the DosPlotter
    plotter = DosPlotter()

    # Add the element projected DOS to the DOSplotter
    for specie in species:
        plotter.add_dos(specie.name, element_dos[specie])

    plotter.show()
