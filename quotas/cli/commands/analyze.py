# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os
import matplotlib.pyplot as plt

from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.vasp.outputs import Locpot, Vasprun

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

    locpot = Locpot.from_file(os.path.join(directory, "LOCPOT"))
    vasprun = Vasprun(os.path.join(directory, "vasprun.xml"))

    average_potential = locpot.get_average_along_axis(2)
    fermi_energy = vasprun.efermi

    if plot_potential:
        plt.plot(average_potential)
        pass #TODO Finish this part

    work_function = max(average_potential) - fermi_energy
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
