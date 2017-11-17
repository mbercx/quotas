
import os
import string

from quotas.calculation import slabRelaxSet, slabWorkFunctionSet,\
    find_suitable_kpar, find_irr_kpoints
from quotas.slab import find_atomic_layers

from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.vasp.inputs import Kpoints, Incar

"""
Setup scripts for the calculations of the quotas package.

"""

# Parameters
MAX_KPAR = 30
DFT_FUNCTIONAL = "PBE_54"
#TODO Find way to make potential setting more user friendly

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
    if not "magmom" in bulk_structure.site_properties.keys():

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

        slab_structure = Structure(slab.lattice, slab.species, slab.frac_coords,
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

    # #TODO Naming is difficult because of the lack of a jsonable Slab class. Fix this after fixing that problem.

    relax_dir = os.path.join(
        current_dir,slab_file.strip(
            slab_structure.composition.reduced_formula + "_"
        ).strip(".json"))

    # Write the input files to the calculation directory
    geo_optimization.write_input(relax_dir)

    if verbose:
        print("Written input files to " + relax_dir)


def wf(relax_dir, k_product):
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

    # Add the KPAR tag to the INCAR file
    kpar(directory=calculation_dir, max_kpar=MAX_KPAR, add_kpar=True)

def dos(relax_dir, k_product):
    """
    Set up the work function calculation based on the output of the geometry
    optimization.

    """
    relax_dir = os.path.abspath(relax_dir)

    # Set up the calculation
    DOS_calc = slabWorkFunctionSet.from_relax_calc(relax_dir=relax_dir,
                                                   k_product=k_product)

    calculation_dir = os.path.join(os.path.split(relax_dir)[0], "dos")

    # Write the input files of the calculation
    DOS_calc.write_input(calculation_dir)

    # Add the KPAR tag to the INCAR file
    kpar(directory=calculation_dir, max_kpar=MAX_KPAR, add_kpar=True)


def slab_setup(bulk_file, miller_indices, thickness, vacuum, fix_part,
               fix_thickness, verbose):
    """
    Set up all the calculations to study the slab of a structure.

    Args:
        bulk_file (string): Name of the structure file.
        miller_indices (list): The miller indices of the surface.
        thickness (float):
        vacuum (float):
        fix (string):

    Returns:

    """

    if verbose:
        print("Importing bulk structure...")

    # Import the bulk structure
    bulk_structure = Structure.from_file(bulk_file)

    # If no magnetic configuration is given, start the calculation in a
    # non-magnetic state.
    if not "magmom" in bulk_structure.site_properties.keys():

        if verbose:
            print("No magnetic configuration found. Adding magmom = 0 for all "
                  "sites.")

        bulk_structure.add_site_property("magmom",
                                         [0]*len(bulk_structure.sites))

    #TODO write checks for oxidation states etc

    if verbose:
        print("Generating slab terminations...")

    # Generate the various terminations of the surface
    slabs = SlabGenerator(bulk_structure, miller_indices,
                          thickness, vacuum).get_slabs()

    if verbose:
        print("Number of slabs found = " + str(len(slabs)))
        print("Removing polar slab terminations...")

    # Remove the polar slabs
    for slab in slabs.copy():
        if slab.is_polar():
            slabs.remove(slab)

    if verbose:
        print("Found " + str(len(slabs)) + " non-polar slab termination(s).")

    current_dir = os.path.dirname(".")

    slab_letter_counter = 0

    if verbose:
        print("Setting up geometry optimizations...")

    for slab in slabs:

        n_atomic_layers = len(find_atomic_layers(slab))

        if verbose:
            print("Number of layers in slab: " + str(n_atomic_layers))

        slab.sort(key=lambda site: site.properties["magmom"])
        slab.sort()

        geo_optimization = slabRelaxSet(slab, potcar_functional=DFT_FUNCTIONAL)
        geo_optimization.fix_slab_bulk(thickness=fix_thickness,
                                       part=fix_part)

        slab_letter = string.ascii_lowercase[slab_letter_counter]
        slab_letter_counter += 1

        geo_dir = "".join([str(number) for number in miller_indices]) + "_" \
                    + slab_letter + "_" + str(n_atomic_layers) + "l" \
                    + str(int(vacuum)) + "v"

        relax_dir = "".join([fix_part, "_relax"])

        calculation_dir = os.path.join(current_dir, geo_dir, relax_dir)

        geo_optimization.write_input(calculation_dir)

        # Add the KPAR tag to the INCAR file
        kpar(directory=calculation_dir, max_kpar=MAX_KPAR, add_kpar=True)

        if verbose:
            print("Written input files to " + calculation_dir)


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