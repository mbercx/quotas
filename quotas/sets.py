# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os
import quotas.core as slab

from monty.serialization import loadfn

from pymatgen.core.structure import Structure
from pymatgen.io.vasp.sets import DictSet
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar

"""
Module that defines the various calculations required for the quotas script.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.2"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Apr 2018"

MODULE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "set_configs")
DFT_FUNCTIONAL = "PBE_54"


def _load_yaml_config(filename):
    config = loadfn(os.path.join(MODULE_DIR, "%s.yaml" % filename))
    return config


class bulkRelaxSet(DictSet):
    """
    VASP input set for the bulk relaxation.

    I prefer having my own yaml files to modify, instead of relying on those
    of the pymatgen package and adding incar_settings.

    """
    CONFIG = _load_yaml_config("bulkRelaxSet")

    def __init__(self, structure, **kwargs):
        super(bulkRelaxSet, self).__init__(
            structure, bulkRelaxSet.CONFIG, **kwargs)
        self.kwargs = kwargs


class bulkSCFSet(DictSet):
    CONFIG = _load_yaml_config("bulkSCFSet")

    def __init__(self, structure, k_product, **kwargs):
        super(bulkSCFSet, self).__init__(
            structure, bulkSCFSet.CONFIG, **kwargs)
        self.k_product = k_product
        self.kwargs = kwargs

    @property
    def kpoints(self):
        """
        Set up the k-points for the bulk SCF.

        The number of k-point divisions in each direction is determined by
        the length of the corresponding lattice vector, by making sure the
        product of the number of k-points and the length of the
        corresponding lattice vector is equal to k_product, defined in the
        initialization of the SCF calculation.

        This means that in order to get a k-
        point spacing of ~0.05
        angstrom^{-1} along the reciprocal lattice vectors, you need to put
        the k_product equal to 20.

        Returns:
            :class: pymatgen.io.vasp.inputs.Kpoints

        """
        # Use k_product to set up kpoints array
        abc = self.structure.lattice.abc
        kpt_calc = [int(self.k_product / abc[0] + 0.5),
                    int(self.k_product / abc[1] + 0.5),
                    int(self.k_product / abc[2] + 0.5)]

        kpoints = Kpoints.gamma_automatic(kpts=kpt_calc)

        return kpoints

    @staticmethod
    def from_relax_calc(relax_dir, k_product, **kwargs):
        """
        Set up the SCF calculation based on the output of the geometry
        optimization.

        """
        relax_dir = os.path.abspath(relax_dir)

        # TODO this can be made more general
        # Obtain the structure from the CONTCAR file of the VASP calculation
        try:
            structure = Structure.from_file(os.path.join(relax_dir, "CONTCAR"))
        except FileNotFoundError:
            structure = Structure.from_file(os.path.join(relax_dir,
                                                         "CONTCAR.vasp"))

        # Try to initialize the magnetic configuration in the same way as for
        # the geometry optimization
        incar = Incar.from_file(os.path.join(relax_dir, "INCAR"))
        try:
            magmom = incar["MAGMOM"]
        except KeyError:
            # If no magnetic moment is present, set it to zero
            magmom = [0] * len(structure.sites)

        structure.add_site_property("magmom", magmom)

        return bulkSCFSet(structure=structure,
                          k_product=k_product,
                          potcar_functional=DFT_FUNCTIONAL,
                          **kwargs)


class slabRelaxSet(DictSet):
    """
    A VASP input set that is used to optimize a slab structure.

    Reason for implementing: I did not immediately find a VASP input set in
    pymatgen that allowed for selective dynamics of the slab.

    """

    CONFIG = _load_yaml_config("slabRelaxSet")

    def __init__(self, structure, k_product=50, **kwargs):
        super(slabRelaxSet, self).__init__(structure=structure,
                                           config_dict=slabRelaxSet.CONFIG,
                                           **kwargs)
        self.k_product = k_product
        self.selective_dynamics = None
        self.kwargs = kwargs

    def fix_slab_bulk(self, thickness, method="layers", part="center"):
        """
        Fix atoms of the slab to represent the bulk of the material. Which atoms are
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
        """

        self.selective_dynamics = slab.fix_slab_bulk(self.poscar, thickness,
                                                     method, part)

    @property
    def poscar(self):
        """
        Similar to the standard POSCAR, but with selective dynamics added.

        Returns:
            :class: pymatgen.io.vasp.inputs.Poscar
        """
        return Poscar(self.structure,
                      selective_dynamics=self.selective_dynamics)

    @property
    def kpoints(self):
        """
        Sets up the k-points for the slab relaxation.

        For slabs, the number of
        k-points corresponding to the third reciprocal lattice vector is set
        to 1. The number of k-point divisions in the other two directions is
        determined by the length of the corresponding lattice vector, by making
        sure the product of the number of k-points and the length of the
        corresponding lattice vector is equal to k_product, defined in the
        initialization of the slab calculation.

        Returns:
            :class: pymatgen.io.vasp.inputs.Kpoints

        """

        # Use k_product to calculate kpoints
        abc = self.structure.lattice.abc
        kpt_calc = [int(self.k_product / abc[0] + 0.5),
                    int(self.k_product / abc[1] + 0.5), 1]

        kpoints = Kpoints.gamma_automatic(kpts=kpt_calc)

        return kpoints


class slabWorkFunctionSet(DictSet):
    """
    A VASP input set that can be used to calculate the work function of a slab.
    """

    CONFIG = _load_yaml_config("slabWorkFunctionSet")

    def __init__(self, structure, k_product=50, **kwargs):
        super(slabWorkFunctionSet,
              self).__init__(structure=structure,
                             config_dict=slabWorkFunctionSet.CONFIG,
                             **kwargs)
        self.k_product = k_product
        self.kwargs = kwargs

    @property
    def kpoints(self):
        """
        Sets up the k-points for the slab relaxation.

        For slabs, the number of
        k-points corresponding to the third reciprocal lattice vector is set
        to 1. The number of k-point divisions in the other two directions is
        determined by the length of the corresponding lattice vector, by making
        sure the product of the number of k-points and the length of the
        corresponding lattice vector is equal to k_product, defined in the
        initialization of the slab calculation.

        Returns:
            :class: pymatgen.io.vasp.inputs.Kpoints

        """

        # Use k_product to calculate kpoints
        abc = self.structure.lattice.abc
        kpt_calc = [int(self.k_product / abc[0] + 0.5),
                    int(self.k_product / abc[1] + 0.5), 1]

        kpoints = Kpoints.gamma_automatic(kpts=kpt_calc)

        return kpoints

    @staticmethod
    def from_relax_calc(relax_dir, k_product, **kwargs):
        """
        Set up the calculation based on the output of the geometry
        optimization.

        """
        relax_dir = os.path.abspath(relax_dir)

        # TODO this can be made more general
        # Obtain the structure from the CONTCAR file of the VASP calculation
        try:
            structure = Structure.from_file(os.path.join(relax_dir, "CONTCAR"))
        except FileNotFoundError:
            structure = Structure.from_file(os.path.join(relax_dir,
                                                         "CONTCAR.vasp"))

        # Initialize the magnetic configuration in the same way as for the
        # geometry optimization
        incar = Incar.from_file(os.path.join(relax_dir, "INCAR"))
        magmom = incar["MAGMOM"]
        structure.add_site_property("magmom", magmom)

        return slabWorkFunctionSet(structure=structure,
                                   k_product=k_product,
                                   potcar_functional=DFT_FUNCTIONAL,
                                   **kwargs)


class slabWorkFunctionHSESet(DictSet):
    """
    A VASP input set that can be used to calculate the work function of a slab.
    """

    CONFIG = _load_yaml_config("slabWorkFunctionHSESet")

    def __init__(self, structure, k_product=50, **kwargs):
        super(slabWorkFunctionHSESet,
              self).__init__(structure=structure,
                             config_dict=slabWorkFunctionHSESet.CONFIG,
                             **kwargs)
        self.k_product = k_product
        self.kwargs = kwargs

    @property
    def kpoints(self):
        """
        Sets up the k-points for the slab relaxation.

        For slabs, the number of
        k-points corresponding to the third reciprocal lattice vector is set
        to 1. The number of k-point divisions in the other two directions is
        determined by the length of the corresponding lattice vector, by making
        sure the product of the number of k-points and the length of the
        corresponding lattice vector is equal to k_product, defined in the
        initialization of the slab calculation.

        Returns:
            :class: pymatgen.io.vasp.inputs.Kpoints

        """

        # Use k_product to calculate kpoints
        abc = self.structure.lattice.abc
        kpt_calc = [int(self.k_product / abc[0] + 0.5),
                    int(self.k_product / abc[1] + 0.5), 1]

        kpoints = Kpoints.gamma_automatic(kpts=kpt_calc)

        return kpoints

    @staticmethod
    def from_relax_calc(relax_dir, k_product, **kwargs):
        """
        Set up the calculation based on the output of the geometry
        optimization.

        """
        relax_dir = os.path.abspath(relax_dir)

        # TODO this can be made more general
        # Obtain the structure from the CONTCAR file of the VASP calculation
        try:
            structure = Structure.from_file(os.path.join(relax_dir, "CONTCAR"))
        except FileNotFoundError:
            structure = Structure.from_file(os.path.join(relax_dir,
                                                         "CONTCAR.vasp"))

        # Initialize the magnetic configuration in the same way as for the
        # geometry optimization
        incar = Incar.from_file(os.path.join(relax_dir, "INCAR"))
        magmom = incar["MAGMOM"]
        structure.add_site_property("magmom", magmom)

        return slabWorkFunctionHSESet(structure=structure,
                                      k_product=k_product,
                                      potcar_functional=DFT_FUNCTIONAL,
                                      **kwargs)
