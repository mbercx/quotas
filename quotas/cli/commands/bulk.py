import os
import subprocess

from monty.serialization import loadfn
from quotas.sets import bulkRelaxSet, bulkSCFSet

from pymatgen.core import Structure
from pymatgen.io.vasp.outputs import Vasprun

"""
Module that defines the commands for setting up bulk calculations.
"""

MODULE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "../../set_configs")
DFT_FUNCTIONAL = "PBE_54"


def _load_yaml_config(filename):
    config = loadfn(os.path.join(MODULE_DIR, "%s.yaml" % filename))
    return config


def relax(bulk_file, is_metal=False, hse_calc=False, verbose=False):

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

    if hse_calc:

        hse_config = _load_yaml_config("HSESet")
        user_incar_settings.update(hse_config["INCAR"])

        # Set up the calculation directory
        current_dir = os.path.dirname(".")
        calculation_dir = os.path.join(current_dir, "bulk", "hse_relax")

    else:

        dftu_config = _load_yaml_config("DFTUSet")
        user_incar_settings.update(dftu_config["INCAR"])

        # Set up the calculation directory
        current_dir = os.path.dirname(".")
        calculation_dir = os.path.join(current_dir, "bulk", "dftu_relax")

    # For metals, add some Methfessel Paxton smearing
    if is_metal:
        user_incar_settings.update({"ISMEAR": 1, "SIGMA": 0.2})


    # Set up the geometry optimization
    geo_optimization = bulkRelaxSet(structure=bulk_structure,
                                    user_incar_settings=user_incar_settings,
                                    potcar_functional=DFT_FUNCTIONAL)

    # Write the input files to the geo optimization directory
    geo_optimization.write_input(calculation_dir)

    if verbose:
        print("Written input to " + calculation_dir)


def dos(relax_dir, k_product, hse_calc=False):
    """
    Set up the work function calculation based on the output of the geometry
    optimization.

    """
    relax_dir = os.path.abspath(relax_dir)

    relax_out = Vasprun(os.path.join(relax_dir, "vasprun.xml"))

    # Make sure enough bands are included in the DOS calculation for the
    # quotas model.
    nbands = relax_out.parameters["NBANDS"]
    while 2*relax_out.parameters['NELECT'] + 4*len(
            relax_out.initial_structure.sites) > nbands:
        nbands += relax_out.parameters["NBANDS"]

    # Add settings for the DOS calculation
    user_incar_settings = {"NEDOS": 2000, "NBANDS":nbands, "LORBIT":11}

    if hse_calc:

        hse_config = _load_yaml_config("HSESet")
        user_incar_settings.update(hse_config["INCAR"])

        # Set up the calculation directory
        calculation_dir = os.path.join(os.path.split(relax_dir)[0], "hse_dos")

    else:

        # TODO check if using the CHGCAR of scf actually works

        dftu_config = _load_yaml_config("DFTUSet")
        user_incar_settings.update(dftu_config["INCAR"])

        # Set up the calculation directory
        calculation_dir = os.path.join(os.path.split(relax_dir)[0], "dftu_dos")

    # Set up the calculation
    dos_calculation = bulkSCFSet.from_relax_calc(
        relax_dir=relax_dir,
        k_product=k_product,
        user_incar_settings=user_incar_settings
    )

    # Write the input files of the calculation
    dos_calculation.write_input(calculation_dir)

    #if not hse_calc:

        # Copy the charge density from the geometry optimization
        #subprocess.call(["cp", os.path.join(relax_dir, "CHGCAR"),
                         #os.path.join(calculation_dir)])


def optics(relax_dir, k_product, hse_calc, is_metal, verbose):

    relax_dir = os.path.abspath(relax_dir)

    relax_out = Vasprun(os.path.join(relax_dir, "vasprun.xml"))

    nbands = relax_out.parameters["NBANDS"] * 3

    # Add the INCAR settings for the dielectric function calculation
    user_incar_settings = {"LOPTICS":True, "NEDOS":2000, "NBANDS": nbands}

    if hse_calc:

        hse_config = _load_yaml_config("HSESet")
        user_incar_settings.update(hse_config["INCAR"])

        # Set up the calculation directory
        calculation_dir = os.path.join(os.path.split(relax_dir)[0],
                                       "hse_optics")

    else:

        dftu_config = _load_yaml_config("DFTUSet")
        user_incar_settings.update(dftu_config["INCAR"])

        # Set up the calculation directory
        calculation_dir = os.path.join(os.path.split(relax_dir)[0],
                                       "dftu_optics")

    # For metals, add some Methfessel Paxton smearing
    if is_metal:

        if verbose:
            print("Metal option detected. Adding good amount of gaussian "
                  "smearing.")

        user_incar_settings.update({"ISMEAR": 0, "SIGMA": 0.3})

    # Set up the calculation
    if verbose:
        print("Setting up calculation...")

    diel_calculation = bulkSCFSet.from_relax_calc(
        relax_dir=relax_dir,
        k_product=k_product,
        user_incar_settings=user_incar_settings
    )

    # Write the input files to the calculation directory
    diel_calculation.write_input(calculation_dir)

    if verbose:
        print("Written input to " + calculation_dir)
