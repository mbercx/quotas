
import os
import quotas.slab as slab

from monty.serialization import loadfn

from pymatgen.io.vasp.sets import DictSet
from pymatgen.io.vasp.inputs import Poscar

"""
Package that defines the various calculations required for the quotas script.

"""

MODULE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "set_configs")

def _load_yaml_config(fname):
    config = loadfn(os.path.join(MODULE_DIR, "%s.yaml" % fname))
    # config["INCAR"].update(loadfn(os.path.join(MODULE_DIR,
    #                                            "VASPIncarBase.yaml")))
    return config

class slabRelaxSet(DictSet):
    """
    A VASP input set that is used to optimize a slab structure.

    Reason for implementing: I did not immediately find a VASP input set in
    pymatgen that allowed for selective dynamics of the slab.

    """

    CONFIG = _load_yaml_config("slabRelaxSet")

    def __init__(self, structure, **kwargs):
        super(slabRelaxSet, self).__init__(structure, slabRelaxSet.CONFIG,
                                           **kwargs)
        self.kwargs = kwargs
        self.selective_dynamics = None

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

        Returns:

        """
        return Poscar(self.structure,
                      selective_dynamics=self.selective_dynamics)



