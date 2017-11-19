# Encoding: UTF-8

import os
import shutil
import string

from quotas.calculation import slabRelaxSet, slabWorkFunctionSet
from quotas.core import find_atomic_layers

from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.vasp.outputs import Vasprun

"""
Setup scripts for the calculations of the quotas package.

"""

# Parameters
MAX_KPAR = 30
DFT_FUNCTIONAL = "PBE_54"
# TODO Find way to make potential setting more user friendly


def setup(bulk_file, miller_indices, thickness, vacuum, verbose):

    if verbose:
        print("Importing bulk structure...")

    # Import the bulk structure
    bulk_structure = Structure.from_file(bulk_file)

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
                                         [0]*len(bulk_structure.sites))

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
        slab_structure.to(fmt="json", filename=slab_file+".json")


def relax(slab_file, fix_part, fix_thickness, verbose):

    if verbose:
        print("Reading structure from file...")

    slab_structure = Structure.from_file(slab_file)

    if verbose:
        print("Setting up calculation...")

    geo_optimization = slabRelaxSet(slab_structure,
                                    potcar_functional=DFT_FUNCTIONAL)
    geo_optimization.fix_slab_bulk(thickness=fix_thickness,
                                   part=fix_part)

    current_dir = os.path.dirname(".")

    # TODO Naming is difficult because of the lack of a jsonable Slab class. Fix this after fixing that problem.

    relax_dir = os.path.join(
        current_dir, slab_file.strip(
            slab_structure.composition.reduced_formula + "_"
        ).strip(".json"))

    # Write the input files to the calculation directory
    geo_optimization.write_input(relax_dir)

    if verbose:
        print("Written input files to " + relax_dir)


def wf(relax_dir, k_product, hse):
    """
    Set up the work function calculation based on the output of the geometry
    optimization.

    """
    relax_dir = os.path.abspath(relax_dir)

    # Set up the calculation
    work_function_calc = \
        slabWorkFunctionSet.from_relax_calc(relax_dir, k_product=k_product)

    calculation_dir = os.path.join(os.path.split(relax_dir)[0], "work_function")

    # Write the input files of the calculation
    work_function_calc.write_input(calculation_dir)


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
        print("Sorry, not implemented yet.")
    else:

        # Use the charge density from the geometry optimization
        dos_incar["ICHARG"] = 11

        # Set up the calculation
        dos_calc = slabWorkFunctionSet.from_relax_calc(
            relax_dir=relax_dir,
            k_product=k_product,
            user_incar_settings=dos_incar
        )

        # Set up the calculation directory
        calculation_dir = os.path.join(os.path.split(relax_dir)[0], "DFTU_dos")

    # Write the input files of the calculation
    dos_calc.write_input(calculation_dir)

    if not hse_calc:

        # Copy the charge density from the geometry optimization
        shutil.copy(os.path.join(relax_dir, "CHGCAR"),
                    os.path.join(calculation_dir))


