
import os
import matplotlib.pyplot as plt

from pymatgen.io.vasp.outputs import Locpot, Vasprun

"""
Module that defines the analysis tools for the quotas package.

"""


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
    print("Work functino = " + work_function)