# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import pdb

import os

from numpy.random import choice

from quotas.core import WorkFunctionData

from quotas.core import find_irr_kpoints, find_suitable_kpar
from pymatgen.core import Structure
from pymatgen.io.vasp.outputs import Vasprun, UnconvergedVASPWarning
from pymatgen.io.vasp.inputs import Incar, Kpoints


"""
Utility commands. Basically something that I can't find a place for elsewhere.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.2"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Apr 2018"


def kpar(directory, max_kpar, add_kpar):
    """

    :return:
    """
    input_dir = os.path.abspath(directory)
    structure = Structure.from_file(os.path.join(input_dir, "POSCAR"))
    kpoints = Kpoints.from_file(os.path.join(input_dir, "KPOINTS"))

    suggested_kpar = str(find_suitable_kpar(structure, kpoints, max_kpar))
    print("Suggested KPAR based on divisors of the number of kpoints = " +
          suggested_kpar)

    if add_kpar:
        print("Adding KPAR tag to INCAR file.")

        try:
            incar = Incar.from_file(os.path.join(directory, "INCAR"))
        except FileNotFoundError:
            raise FileNotFoundError("The INCAR file is not found in the "
                                    "directory.")

        incar["KPAR"] = suggested_kpar
        incar.write_file(os.path.join(directory, "INCAR"))


def nkp(directory):
    """

    Args:
        directory:

    Returns:

    """

    # TODO Seems to fail. Could be related to the fact that the symmetry
    # tolerances are different.

    input_dir = os.path.abspath(directory)
    structure = Structure.from_file(os.path.join(input_dir, "POSCAR"))
    kpoints = Kpoints.from_file(os.path.join(input_dir, "KPOINTS"))

    print("Number of irreducible kpoints = " +
          str(find_irr_kpoints(structure, kpoints)))


def check_run(directory):
    """
    Check a VASP run to see if it has completed successfully.

    Args:
        directory:

    Returns:

    """
    import warnings
    warnings.filterwarnings("error")

    # TODO This is pretty slow because the whole vasprun has to be parsed.
    # Improve function somehow.

    try:
        out = Vasprun(os.path.join(directory, "vasprun.xml"))
    except FileNotFoundError:
        print("No vasprun.xml found in run directory.")
    except UnconvergedVASPWarning:
        print("Calculation has not converged successfully.")


def process_output(directory, calculation="dos"):
    """


    Args:
        directory:

    Returns:

    """

    wf = WorkFunctionData.from_output(
        poscar_filename=os.path.join(directory, "POSCAR"),
        locpot_filename=os.path.join(directory, "LOCPOT"),
        outcar_filename=os.path.join(directory, "OUTCAR")
    )

    wf.to(os.path.join(directory, "wf_data.json"))

def generate_alloy(root_structure_file, alloy_element, mixing_ratio,
                   number_random):
    """

    Args:
        root_structure:
        alloy_element:
        mixing_ratio:
        number_random:

    Returns:

    """

    alloy = Structure.from_file(root_structure_file)

    random_list = []

    # Make a list of random combinations of site indices
    while len(random_list) < number_random:

        random_indices = list(choice(
            list(range(len(alloy))), round(len(alloy) * mixing_ratio), False
        ))

        if random_indices not in random_list:
            random_list.append(random_indices)

    # Get the current working directory
    directory = os.getcwd()

    # Set up the list of random structures
    for i, index_list in zip(range(len(random_list)), random_list):

        random_structure = alloy.copy()

        for site_index in index_list:
            random_structure.replace(site_index, alloy_element)

        filename = os.path.basename(root_structure_file).strip(
            ".json").strip(".cif") + "_rand" + str(i) + ".cif"

        random_structure.sort()
        random_structure.to(
            fmt="cif", filename=os.path.join(directory, filename)
        )
