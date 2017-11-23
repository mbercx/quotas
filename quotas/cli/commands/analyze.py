
import os
import matplotlib.pyplot as plt

from pymatgen.io.vasp.outputs import Locpot

"""
Module that defines the analysis tools for the quotas package.

"""


def wf(directory, plot_potential=False):

    directory = os.path.abspath(directory)

    locpot_file = os.path.join(directory, "LOCPOT")

    locpot = Locpot.from_file(locpot_file)

    average_potential = locpot.get_average_along_axis(2)

    if plot_potential:
        plt.plot(average_potential)
        pass #TODO Finish this part

    return max(average_potential)