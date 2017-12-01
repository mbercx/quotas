import os
import subprocess

from quotas.calculation import bulkRelaxSet, bulkDosSet, bulkDosHSESet

from pymatgen.core import Structure
from pymatgen.io.vasp.outputs import Vasprun

"""
Module that defines the commands for setting up bulk calculations.
"""

DFT_FUNCTIONAL = "PBE_54"


def relax(bulk_file, is_metal, verbose):

    if verbose:
        print("Reading structure from file...")

    # Read the bulk structure
    bulk_structure = Structure.from_file(bulk_file)

    # If no magnetic configuration is given, start the calculation in a
    # non-magnetic state.
    if "magmom" not in bulk_structure.site_properties.keys():

        if verbose:
            print("No magnetic configuration found. Adding magmom = 0 for all "
                  "sites.")

        bulk_structure.add_site_property("magmom",
                                         [0] * len(bulk_structure.sites))

    if verbose:
        print("Setting up calculation...")

    user_incar_settings = {}

    # For metals, add some Methfessel Paxton smearing
    if is_metal:
        user_incar_settings = {"ISMEAR": 1, "SIGMA": 0.2}


    # Set up the geometry optimization
    geo_optimization = bulkRelaxSet(structure=bulk_structure,
                                    user_incar_settings=user_incar_settings,
                                    potcar_functional=DFT_FUNCTIONAL)

    # Set up the geometry optimization directory
    current_dir = os.path.dirname(".")

    relax_dir = os.path.join(current_dir, "bulk", "relax")

    # Write the input files to the geo optimization directory
    geo_optimization.write_input(relax_dir)

    if verbose:
        print("Written input to " + relax_dir)


def dos(relax_dir, k_product, hse_calc=False):
    """
    Set up the work function calculation based on the output of the geometry
    optimization.

    """
    relax_dir = os.path.abspath(relax_dir)

    relax_out = Vasprun(os.path.join(relax_dir, "vasprun.xml"))

    nbands = relax_out.parameters["NBANDS"]*3

    # Add some typical extra settings for the DOS calculation
    dos_incar = {"NEDOS": 2000, "NBANDS":nbands}

    if hse_calc:

        # Set up the calculation
        dos_calc = bulkDosHSESet.from_relax_calc(
            relax_dir=relax_dir,
            k_product=k_product,
            user_incar_settings=dos_incar
        )

        # Set up the calculation directory
        calculation_dir = os.path.join(os.path.split(relax_dir)[0], "hse_dos")

    else:

        # Use the charge density from the geometry optimization
        dos_incar["ICHARG"] = 11

        # Set up the calculation
        dos_calc = bulkDosSet.from_relax_calc(
            relax_dir=relax_dir,
            k_product=k_product,
            user_incar_settings=dos_incar
        )

        # Set up the calculation directory
        calculation_dir = os.path.join(os.path.split(relax_dir)[0], "dftu_dos")

    # Write the input files of the calculation
    dos_calc.write_input(calculation_dir)

    if not hse_calc:

        # Copy the charge density from the geometry optimization
        subprocess.call(["cp", os.path.join(relax_dir, "CHGCAR"),
                         os.path.join(calculation_dir)])

def diel(relax_dir, is_metal, hse_calc, verbose):

    relax_dir = os.path.abspath(relax_dir)

    relax_out = Vasprun(os.path.join(relax_dir, "vasprun.xml"))

    nbands = relax_out.parameters["NBANDS"] * 3

    if hse_calc:

        # Set up the calculation
        dos_calc = bulkDosHSESet.from_relax_calc(
            relax_dir=relax_dir,
            k_product=k_product,
            user_incar_settings=dos_incar
        )

        # Set up the calculation directory
        calculation_dir = os.path.join(os.path.split(relax_dir)[0], "hse_dos")

    else:

        # Set up the calculation
        dos_calc = bulkDosSet.from_relax_calc(
            relax_dir=relax_dir,
            k_product=k_product,
            user_incar_settings=dos_incar
        )

        # Set up the calculation directory
        calculation_dir = os.path.join(os.path.split(relax_dir)[0], "dftu_dos")



    if verbose:
        print("Setting up calculation...")

    user_incar_settings = {"NEDOS": 2000, "NBANDS": nbands}

    # For metals, add some Methfessel Paxton smearing
    if is_metal:
        user_incar_settings.update({"ISMEAR": 1, "SIGMA": 0.3})


    # Set up the geometry optimization
    geo_optimization = bulkRelaxSet(structure=bulk_structure,
                                    user_incar_settings=user_incar_settings,
                                    potcar_functional=DFT_FUNCTIONAL)

    # Set up the geometry optimization directory
    current_dir = os.path.dirname(".")

    relax_dir = os.path.join(current_dir, "bulk", "diel")

    # Write the input files to the geo optimization directory
    geo_optimization.write_input(relax_dir)

    if verbose:
        print("Written input to " + relax_dir)
