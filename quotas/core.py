# Encoding: UTF-8

from pymatgen.core.surface import SlabGenerator
import string


"""
A set of methods to aid in the setup of slab calculations. 

"""

def fix_slab_bulk(poscar, thickness, method="layers", part="center"):
    """
    Fix atoms of a slab to represent the bulk of the material. Which atoms are
    fixed depends on whether the user wants to fix one side or the center, and
    how exactly the part of the slab is defined.

    Args:
        poscar (:class: `pymatgen.io.vasp,outputs.Poscar`): The poscar file of the slab

        thickness (float): The thickness of the fixed part of the slab,
            expressed in number of layers or Angstroms, depending on the
            method.

        method (string): How to define the thickness of the part of the slab
            that is fixed:

                "layers" (default): Fix a set amount of layers. The layers are
                found using the 'find_atomic_layers' method.
                "angstroms": Fix a part of the slab of a thickness defined in
                angstroms.

        part (string): Which part of the slab to fix:

                "center" (default): Fix the atoms at the center of the slab.


    Returns:
        selective dynamics (Nx3 array): bool values for selective dynamics,
            where N is number of sites. Defaults to None.

    """
    if method == "layers":

        atomic_layers = find_atomic_layers(poscar.structure)

        if part == "center":

            # Even number of layers
            if len(atomic_layers)%2 == 0:

                # Check if the user requested an odd number of layers for the
                # fixed part of the slab
                if thickness%2 == 1:

                    print("Found an even number of layers, but the user " +
                          "requested an odd number of fixed layers. Adding "
                          "one layer to the fixed part of the slab.")
                    thickness += 1

            # Odd number of layers
            if len(atomic_layers)%2 == 1:

                # Check if the user requested an even number of layers for the
                # fixed part of the slab
                if thickness % 2 == 0:
                    print("Found an odd number of layers, but the user " +
                          "requested an even number of fixed layers. Adding "
                          "one layer to the fixed part of the slab.")
                    thickness += 1

            # Calculate the number of layers to optimiz on each site
            n_optimize_layers = int((len(atomic_layers) - thickness)/2)

            # Take the fixed layers from the atomic layers of the slab
            fixed_layers = atomic_layers[n_optimize_layers: n_optimize_layers +
                                                            thickness]

        else:
            raise NotImplementedError("Requested part is not implemented " +
                                      "(yet).")
            #TODO Implement oneside

        fixed_sites = [site for layer in fixed_layers for site in layer]

        # Set up the selective dynamics property

        selective_dynamics = []

        for site in poscar.structure.sites:
            if site in fixed_sites:
                selective_dynamics.append([False, False, False])
            else:
                selective_dynamics.append([True, True, True])

        return selective_dynamics

    else:
        raise NotImplementedError("Requested method is not implemented (yet).")
        #TODO Implement angstrom

def find_atomic_layers(structure, layer_tol=1e-2):
    """
    Determines the atomic layers in the c-direction of a structure. Usually
    used to calculate the number of atomic layers in a slab. Note that as long
    as a site is "close enough" to ONE other site of a layer (determined by the
    'layer_tol' variable), it will be added to that layer. Another option would
    be to demand that the distance is smaller than 'layer_tol' for ALL sites of
    the layer, but then the division in layers could depend on the order of the
    sites.

    Args:
        structure (pymatgen.core.structure.Structure): Structure for which to
            analyze the layers.
        layer_tol (float): Tolerance for the difference between the third
            fractional coordinate of two sites in the same layer.

    Returns:
        (list) List of the atomic layers, sorted by their position in the c-direction.
        Each atomic layer is also represented by a list of sites.
    """

    atomic_layers = []

    for site in structure.sites:

        is_in_layer = False

        # Check to see if the site is in a layer that is already in our list
        for layer in atomic_layers:

            # Compare the third fractional coordinate of the site with that of
            # the atoms in the considered layer
            for atom_site in layer.copy():
                 if abs(atom_site.frac_coords[2]
                                - site.frac_coords[2]) < layer_tol:
                     is_in_layer = True
                     layer.append(site)
                     break # Break out of the loop, else the site is added
                           # multiple times

        # If the site is not found in any of the atomic layers, create a new
        # atomic layer
        if is_in_layer == False:
            atomic_layers.append([site,])

    # Sort the atomic layers
    atomic_layers.sort(key=lambda layer: layer[0].frac_coords[2])

    return atomic_layers

def write_all_slab_terminations(structure, miller_indices, min_slab_size,
                                min_vacuum_size):
    """
    Writes the POSCAR files for all slab terminations.

    Args:
        structure:
        miller_indices:
        min_slab_size:
        min_vacuum_size:

    Returns:
        None
    """

    slab_gen = SlabGenerator(structure, miller_indices, min_slab_size,
                             min_vacuum_size)
    slabs = slab_gen.get_slabs()

    letter_counter = 0

    for slab in slabs:

        slab.sort(reverse=True)

        slab_letter = string.ascii_lowercase[letter_counter]
        letter_counter += 1

        filename = "".join([str(number) for number in miller_indices]) + "_" \
                    + slab_letter + "_" + str(min_slab_size) + "l_POSCAR.vasp"

        slab.to(fmt='vasp', filename=filename)

        print("Written slab structure to " + filename)
        print("Slab is symmetric = " + str(slab.is_symmetric()))
        print("Slab is polar = " + str(slab.is_polar()))


def find_suitable_kpar(structure, kpoints, max_kpar=30):
    """

    :param structure:
    :param kpoints:
    :return:
    """

    spg = SpacegroupAnalyzer(structure)

    kpar = len(spg.get_ir_reciprocal_mesh(kpoints.kpts))
    divisors = generate_divisors(kpar)

    while kpar > max_kpar:
        kpar = next(divisors)

    return kpar


def find_irr_kpoints(structure, kpoints):
    """

    :param structure:
    :param kpoints:
    :return:
    """

    spg = SpacegroupAnalyzer(structure)

    #TODO Find and fix bug!
    return len(spg.get_ir_reciprocal_mesh(kpoints.kpts))


def generate_divisors(number):
    """

    Args:
        number:

    Returns:

    """
    divisor = int(number/2)

    while divisor != 0:

        while number%divisor != 0:
            divisor -= 1

        yield divisor
        divisor -= 1
