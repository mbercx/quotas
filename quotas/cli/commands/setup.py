
import os
import string

from quotas.calculation import slabRelaxSet
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator

"""
Setup scripts for the calculations of the quotas package.

"""

def slab_setup(filename, miller_indices, thickness, vacuum, fix_part,
               fix_thickness, verbose):
    """
    Set up all the calculations to study the slab of a structure.

    Args:
        filename (string): Name of the structure file.
        miller_indices (list): The miller indices of the surface.
        thickness (float):
        vacuum (float):
        fix (string):

    Returns:

    """

    if verbose:
        print("Importing bulk structure...")

    # Import the bulk structure
    bulk_structure = Structure.from_file(filename)

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
        print("Found " + str(len(slabs)) + " non-polar slab terminations.")

    current_dir = os.path.dirname(".")

    slab_letter_counter = 0

    for slab in slabs:

        slab.sort(key=lambda site: site.properties["magmom"])
        slab.sort()

        geo_optimization = slabRelaxSet(slab, potcar_functional="PBE_54")
        geo_optimization.fix_slab_bulk(thickness=fix_thickness,
                                       part=fix_part)

        slab_letter = string.ascii_lowercase[slab_letter_counter]
        slab_letter_counter += 1

        geo_dir = "".join([str(number) for number in miller_indices]) + "_" \
                    + slab_letter + "_" + str(thickness) + "l"

        relax_dir = "".join([fix_part, "_relax"])

        geo_optimization.write_input(os.path.join(current_dir, geo_dir,
                                                  relax_dir))






