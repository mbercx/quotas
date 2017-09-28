
import os
import string

from quotas.calculation import slabRelaxSet
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator

"""
Setup scripts for the calculations of the quotas package.

"""

def slab_setup(filename, miller_indices, thickness, vacuum, fix_part,
               fix_thickness):
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
    # Import the bulk structure
    bulk_structure = Structure.from_file(filename)

    #TODO write checks for oxidation states etc

    # Generate the various terminations of the surface
    slabs = SlabGenerator(bulk_structure, miller_indices,
                          thickness, vacuum).get_slabs()

    # Remove the polar slabs
    for slab in slabs.copy():
        if slab.is_polar():
            slabs.remove(slab)

    current_dir = os.path.dirname(".")

    slab_letter_counter = 0

    for slab in slabs:

        slab.sort(key=lambda site: site.properties["magmom"])
        slab.sort()

        print(slab.site_properties)

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






