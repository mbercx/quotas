# Encoding: UTF-8

import os
import string

import pdb

from quotas.sets import slabRelaxSet, slabWorkFunctionSet, \
    slabWorkFunctionHSESet
from quotas.core import find_atomic_layers

from monty.serialization import loadfn

from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.vasp.outputs import Vasprun, Outcar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

"""
Setup scripts for the calculations of the quotas package.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.2"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Apr 2018"

MODULE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "../../set_configs")
DFT_FUNCTIONAL = "PBE_54"


def _load_yaml_config(filename):
    config = loadfn(os.path.join(MODULE_DIR, "%s.yaml" % filename))
    return config


# Parameters
MAX_KPAR = 30


# TODO Find way to make potential setting more user friendly


def setup(bulk_file, miller_indices, thickness, vacuum, write_cif, verbose):
    """
    Set up the slab from the bulk structure file. Automatically converts the
    bulk structure to the conventional lattice.

    Args:
        bulk_file:
        miller_indices:
        thickness:
        vacuum:
        write_cif:
        verbose:

    Returns:

    """

    if verbose:
        print("Importing bulk structure...")

    # Import the bulk structure
    bulk_structure = Structure.from_file(bulk_file)

    # Convert the structure to the conventional lattice
    spg = SpacegroupAnalyzer(bulk_structure)
    bulk_structure = spg.get_conventional_standard_structure()

    # Check if the provided structure is oxidation state decorated.
    try:
        test = bulk_structure.sites[0].specie.oxi_state
    except AttributeError:
        print("WARNING: No oxidation state decoration found. Cannot properly "
              "determine slab polarity.")

    # If no magnetic configuration is given, start the calculation in a
    # non-magnetic state.
    if "magmom" not in bulk_structure.site_properties.keys():

        if verbose:
            print("No magnetic configuration found. Adding magmom = 0 for all "
                  "sites.")

        bulk_structure.add_site_property("magmom",
                                         [0] * len(bulk_structure.sites))

    if verbose:
        print("Generating slab terminations...")

    # Generate the various terminations of the surface
    slabs = SlabGenerator(bulk_structure, miller_indices,
                          thickness, vacuum).get_slabs()

    if verbose:
        print("Number of slabs found = " + str(len(slabs)))

    # Write the structure files (.json format for site properties)

    if verbose:
        print("Writing structure files...")

    # Keep track of letters for terminations
    slab_letter_counter = 0

    for slab in slabs:

        n_atomic_layers = len(find_atomic_layers(slab))

        if verbose:
            print("Number of layers in slab: " + str(n_atomic_layers))

        slab.sort(key=lambda site: site.properties["magmom"])
        slab.sort()

        slab_letter = string.ascii_lowercase[slab_letter_counter]
        slab_letter_counter += 1

        slab_file = bulk_structure.composition.reduced_formula + "_" \
                    + "".join([str(number) for number in miller_indices]) \
                    + "_" + slab_letter + "_" + str(n_atomic_layers) + "l" \
                    + str(int(vacuum)) + "v"

        # Add an extra tag to the name in case the slab is polar
        if slab.is_polar():
            slab_file += "_polar"

        # Directly writing a slab to a json doesn't seem to work. So as a
        # workaround, I'll define a structure for each slab with the slab
        # properties.

        # #TODO Change the pymatgen Slab to allow for .json serialization

        # #TODO And allow a slab to be constructed from a .json!

        slab_structure = Structure(slab.lattice, slab.species,
                                   slab.frac_coords,
                                   site_properties=slab.site_properties)
        slab_structure.to(fmt="json", filename=slab_file + ".json")

        if write_cif:
            slab_structure.to(fmt="cif", filename=slab_file + ".cif")


def relax(structure_file, fix_part, fix_thickness, is_metal=False,
          verbose=False):
    """

    Args:
        structure_file:
        fix_part:
        fix_thickness:
        is_metal:
        verbose:

    Returns:
        relax_dir: Full path to the directory where the geometry
        optimization was set up.
    """

    if verbose:
        print("Reading structure from file...")

    slab_structure = Structure.from_file(structure_file)

    # If no magnetic configuration is given, start the calculation in a
    # non-magnetic state.
    if "magmom" not in slab_structure.site_properties.keys():

        if verbose:
            print("No magnetic configuration found. Adding magmom = 0 for all "
                  "sites.")

        slab_structure.add_site_property("magmom",
                                         [0] * len(slab_structure.sites))

    if verbose:
        print("Setting up calculation...")

    if is_metal:

        geo_optimization = slabRelaxSet(slab_structure,
                                        user_incar_settings={"ISMEAR": 1,
                                                             "SIGMA": 0.2},
                                        potcar_functional=DFT_FUNCTIONAL)
        geo_optimization.fix_slab_bulk(thickness=fix_thickness,
                                       part=fix_part)

    else:

        geo_optimization = slabRelaxSet(slab_structure,
                                        potcar_functional=DFT_FUNCTIONAL)
        geo_optimization.fix_slab_bulk(thickness=fix_thickness,
                                       part=fix_part)

    current_dir = os.getcwd()

    # TODO Naming is difficult because of the lack of a jsonable Slab class.
    # Fix this after fixing that problem.

    calculation_dir = os.path.join(
        current_dir, structure_file.strip(
            slab_structure.composition.reduced_formula + "_"
        ).strip(".json").strip(".cif"), fix_part + "_" + "relax")

    if os.path.exists(calculation_dir):
        clean_dir(calculation_dir)

    # Write the input files to the calculation directory
    geo_optimization.write_input(calculation_dir)

    if verbose:
        print("Written input files to " + calculation_dir)

    # Return absolute path to directory of geometry optimization for workflow
    # purposes
    return calculation_dir


def wf(relax_dir, k_product, hse_calc=False):
    """
    Set up the work function calculation based on the output of the geometry
    optimization.

    """
    relax_dir = os.path.abspath(relax_dir)

    user_incar_settings = {}

    if hse_calc:

        hse_config = _load_yaml_config("HSESet")
        user_incar_settings.update(hse_config["INCAR"])

        # Set up the calculation directory
        calculation_dir = os.path.join(os.path.split(relax_dir)[0],
                                       "hse_wf")

    else:

        dftu_config = _load_yaml_config("DFTUSet")
        user_incar_settings.update(dftu_config["INCAR"])

        # Set up the calculation directory
        calculation_dir = os.path.join(os.path.split(relax_dir)[0],
                                       "dftu_wf")

    if os.path.exists(calculation_dir):
        clean_dir(calculation_dir)

    # Set up the calculation
    work_function_calc = slabWorkFunctionSet.from_relax_calc(
        relax_dir=relax_dir,
        k_product=k_product,
        user_incar_settings=user_incar_settings
    )

    # Write the input files of the calculation
    work_function_calc.write_input(calculation_dir)


def dos(relax_dir, k_product, hse_calc=False):
    """
    Set up the work function calculation based on the output of the geometry
    optimization.

    """
    relax_dir = os.path.abspath(relax_dir)

    try:
        relax_vasprun = Vasprun(os.path.join(relax_dir, "vasprun.xml"))

        # Triple the amount of bands compared to the minimum
        nbands = relax_vasprun.parameters["NBANDS"] * 3

    except FileNotFoundError:
        relax_outcar = Outcar(os.path.join(relax_dir, "OUTCAR"))

        relax_outcar.read_pattern(
            {"nbands":  r"\s+k-points\s+NKPTS =\s+[0123456789]+\s+k-points "
                        r"in BZ\s+NKDIM =\s+[0123456789]+\s+number of "
                        r"bands\s+NBANDS=\s+([\.\-\d]+)"},
            postprocess=int)

        # Include a significant number of empty bands
        nbands = relax_outcar.data['nbands'][0][0] * 3

    # Add some typical extra settings for the DOS calculation
    dos_incar = {"NEDOS": 2000, "NBANDS": nbands}

    if hse_calc:

        # Set up the calculation
        dos_calc = slabWorkFunctionHSESet.from_relax_calc(
            relax_dir=relax_dir,
            k_product=k_product,
            user_incar_settings=dos_incar
        )

        # Set up the calculation directory
        calculation_dir = os.path.join(os.path.split(relax_dir)[0], "hse_dos")

    else:

        # Use the charge density from the geometry optimization
        # dos_incar["ICHARG"] = 11

        # Set up the calculation
        dos_calc = slabWorkFunctionSet.from_relax_calc(
            relax_dir=relax_dir,
            k_product=k_product,
            user_incar_settings=dos_incar
        )

        # Set up the calculation directory
        calculation_dir = os.path.join(os.path.split(relax_dir)[0], "dftu_dos")

    if os.path.exists(calculation_dir):
        clean_dir(calculation_dir)

    # Write the input files of the calculation
    dos_calc.write_input(calculation_dir)

    # Return the calculation director for workflow purposes
    return calculation_dir


def clean_dir(directory):
    """
    Clean up a directory so a new calculation can be started in it.

    Right now, it just removes all files from the directory. Might want to
    change that later.

    """
    for file in os.listdir(directory):
        file_path = os.path.join(directory, file)
        try:
            if os.path.isfile(file_path):
                os.remove(file_path)
        except Exception as e:
            print(e)
