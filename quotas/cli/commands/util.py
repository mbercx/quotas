# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os

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

    wf.to(os.path.join("wf_data.json"))





