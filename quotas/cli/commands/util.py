import os

from quotas.core import find_irr_kpoints, find_suitable_kpar
from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Incar, Kpoints


"""
Utility commands. Basically something that I can't find a place for elsewhere.
"""

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

    input_dir = os.path.abspath(directory)
    structure = Structure.from_file(os.path.join(input_dir, "POSCAR"))
    kpoints = Kpoints.from_file(os.path.join(input_dir, "KPOINTS"))

    print("Number of irreducible kpoints = " +
          str(find_irr_kpoints(structure, kpoints)))