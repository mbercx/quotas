# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import os
import subprocess
import shlex

from pymatgen.io.vasp import VaspInput

from custodian import Custodian
from custodian.utils import backup
from custodian.vasp.handlers import VaspErrorHandler, \
    UnconvergedErrorHandler
from custodian.vasp.jobs import VaspJob
from custodian.vasp.interpreter import VaspModder

from fireworks import Firework, LaunchPad, PyTask, FWorker, \
    Workflow
from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter

"""
Workflow setup for the quotas calculations.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.2"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Apr 2018"

# Directory with templates for the Template class
TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            "templates")
VASP_RUN_SCRIPT = "/user/antwerpen/202/vsc20248/local/scripts/job_workflow.sh"
VASP_RUN_COMMAND = "bash /user/antwerpen/202/vsc20248/local/scripts" \
                   "/job_workflow.sh"

# Set up the Launchpad for the workflows
LAUNCHPAD = LaunchPad(host="ds135179.mlab.com", port=35179, name="quotas",
                      username="mbercx", password="quotastests")

def run_vasp(directory):
    """
    Method that simply runs VASP in the directory that is specified. Mainly
    used to set up a PyTask that uses outputs from PyTasks in previous
    FireWorks to run VASP in the appropriate directory.

    Args:
        directory: Absolute path to the directory in which VASP should be run.

    Returns:
        None

    """

    os.chdir(directory)
    subprocess.call(VASP_RUN_SCRIPT)


def run_custodian(directory):
    """
    Run VASP under supervision of a custodian in a certain directory.

    Args:
        directory:

    Returns:

    """

    directory = os.path.abspath(directory)
    os.chdir(directory)

    output = os.path.join(directory, "out")
    vasp_cmd = shlex.split(VASP_RUN_COMMAND)

    # Choose not to use certain error messages to be handled
    error_subset = list(VaspErrorHandler.error_msgs.keys())
    error_subset.remove("brmix")
    vasp_handler = VaspErrorHandler(output_filename=output,
                                    errors_subset_to_catch=error_subset)

    quotas_handler = QuotasErrorHandler(output_filename=output)

    handlers = [quotas_handler,
                UnconvergedErrorHandler(output_filename=output)]

    jobs = [VaspJob(vasp_cmd=vasp_cmd,
                    output_file=output,
                    stderr_file=output)]

    c = Custodian(handlers, jobs, max_errors=3)
    c.run()


def dos_workflow(structure_file, fix_part, fix_thickness, is_metal,
                 k_product, in_custodian):
    """
    Set up a workflow to calculate the DOS of a slab file. Will set up two
    FireWorks:

    "Slab Geometry Optimization": Optimizes the geometry of the slab given
    by the structure file. The user can also specify the part of the slab
    that is fixed using selective dynamics, as well as the number of layers
    to fix.

    "DOS Calculation": Calculates the density of states of the slab, as well as
    the local potential, which can be used to determine the vacuum level.

    After it is set up, the workflow is sent to the LAUNCHPAD.

    Args:
        structure_file (str): Name of the structure file which contains the
            slab geometry.
        fix_part (str): Defines the part of the slab that will remain fixed
            during the geometry optimization.
        fix_thickness (int): Number of atomic layers to fix for the geometry
            optimization.
        is_metal (bool): Specifies whether or not the material is metallic.
            Sets the smearing method to Methfessel-Paxton.
        k_product (int): Determines the density of the k-mesh in the density of
            states calculation. k_product represents the product of the number
        of k-points corresponding to a certain lattice vector with the
            length of that lattice vector.
    """

    # TODO add checks
    # Currently the workflow will be submitted even if the file doesn't exist!

    current_dir = os.getcwd()

    # Set up the geometry optimization from the structure file. All input is
    # provided by the CLI arguments and options. The directory where the
    # geometry optimization is set up is returned and passed as output,
    # so it can be used by Firework children to run the calculation.
    setup_relax = PyTask(func="quotas.cli.commands.slab.relax",
                         kwargs={"structure_file": structure_file,
                                 "fix_part": fix_part,
                                 "fix_thickness": fix_thickness,
                                 "is_metal": is_metal,
                                 "verbose": False},
                         outputs=["relax_dir"]
                         )

    if in_custodian:
        # Run the VASP calculation within a Custodian
        run_relax = PyTask(func="quotas.workflow.run_custodian",
                           inputs=["relax_dir"])
    else:
        # Run the VASP calculation.
        run_relax = PyTask(func="quotas.workflow.run_vasp",
                           inputs=["relax_dir"])

    relax_firework = Firework(tasks=[setup_relax, run_relax],
                              name="Slab Geometry optimization",
                              spec={"_launch_dir": current_dir})

    # Set up the DOS calculation, based on the structure found from the
    # geometry optimization.
    setup_dos = PyTask(func="quotas.cli.commands.slab.dos",
                       inputs=["relax_dir"],
                       kwargs={"k_product": k_product},
                       outputs=["dos_dir"])

    if in_custodian:
        # Run the VASP calculation within a Custodian
        run_dos = PyTask(func="quotas.workflow.run_custodian",
                           inputs=["dos_dir"])
    else:
        # Run the VASP calculation
        run_dos = PyTask(func="quotas.workflow.run_vasp",
                           inputs=["dos_dir"])

    dos_firework = Firework(tasks=[setup_dos, run_dos],
                            name="DOS calculation")

    # Postprocessing
    # TODO Add postprocessing firework

    # Extract the necessary output
    # Calculate the work function -> send to database?

    # Add the workflow to the launchpad
    workflow = Workflow(fireworks=[relax_firework, dos_firework],
                        links_dict={relax_firework: [dos_firework]},
                        name=structure_file + " DOS calculation")

    LAUNCHPAD.add_wf(workflow)


def bulk_optics_workflow(structure_file, is_metal, hse_calc, k_product,
                         in_custodian):
    """
    Sets up a workflow that calculates the dielectric function for the bulk
    structure of a material.

    Returns:

    """

    current_dir = os.getcwd()

    # Set up the geometry optimization from the structure file. All input is
    # provided by the CLI arguments and options. The directory where the
    # geometry optimization is set up is returned and passed as output,
    # so it can be used by Firework children to run the calculation.
    setup_relax = PyTask(func="quotas.cli.commands.bulk.relax",
                         kwargs={"structure_file": structure_file,
                                 "is_metal": is_metal,
                                 "hse_calc": hse_calc,
                                 "verbose": False},
                         outputs=["relax_dir"]
                         )

    if in_custodian:
        # Run the VASP calculation within a Custodian
        run_relax = PyTask(func="quotas.workflow.run_custodian",
                           inputs=["relax_dir"])
    else:
        # Run the VASP calculation
        run_relax = PyTask(func="quotas.workflow.run_vasp",
                           inputs=["relax_dir"])

    relax_firework = Firework(tasks=[setup_relax, run_relax],
                              name="Bulk Geometry optimization",
                              spec={"_launch_dir": current_dir})

    # Set up the dielectric function calculation
    setup_optics = PyTask(func="quotas.cli.commands.bulk.optics",
                          inputs=["relax_dir"],
                          kwargs={"k_product": k_product,
                                  "is_metal": is_metal,
                                  "hse_calc": hse_calc,
                                  "verbose": False},
                          outputs=["optics_dir"])

    if in_custodian:
        # Run the VASP calculation within a Custodian
        run_optics = PyTask(func="quotas.workflow.run_custodian",
                           inputs=["optics_dir"])
    else:
        # Run the VASP calculation
        run_optics = PyTask(func="quotas.workflow.run_vasp",
                           inputs=["optics_dir"])

    optics_firework = Firework(tasks=[setup_optics, run_optics],
                               name="Optics calculation")

    # Add the workflow to the launchpad
    workflow = Workflow(fireworks=[relax_firework, optics_firework],
                        links_dict={relax_firework: [optics_firework]},
                        name=structure_file + " Optics calculation")

    LAUNCHPAD.add_wf(workflow)


def test_custodian(structure_file, fix_part, fix_thickness, is_metal,
                   k_product):
    """
    Testscript for using Custodian to gracefully recover from errors.

    Returns:
        None

    """

    current_dir = os.getcwd()

    # Set up the geometry optimization from the structure file. All input is
    # provided by the CLI arguments and options. The directory where the
    # geometry optimization is set up is returned and passed as output,
    # so it can be used by Firework children to run the calculation.
    setup_relax = PyTask(func="quotas.cli.commands.slab.relax",
                         kwargs={"structure_file": structure_file,
                                 "fix_part": fix_part,
                                 "fix_thickness": fix_thickness,
                                 "is_metal": is_metal,
                                 "verbose": False},
                         outputs=["relax_dir"]
                         )

    # Run the VASP calculation.
    run_relax = PyTask(func="quotas.workflow.run_custodian",
                       inputs=["relax_dir"])

    relax_firework = Firework(tasks=[setup_relax, run_relax],
                              name="Slab Geometry optimization",
                              spec={"_launch_dir": current_dir})

    # -----> Here we would add a check to see if the job completed
    # successfully. If not, we can add another FireWork that makes the
    # necessary adjustments and restarts the calculation.

    # Set up the calculation
    setup_dos = PyTask(func="quotas.cli.commands.slab.dos",
                       inputs=["relax_dir"],
                       kwargs={"k_product": k_product},
                       outputs=["dos_dir"])

    # Run the VASP calculation.
    run_dos = PyTask(func="quotas.workflow.run_custodian",
                     inputs=["dos_dir"])

    dos_firework = Firework(tasks=[setup_dos, run_dos],
                            name="DOS calculation")

    # ----> Here we would add another check...

    ## Firework 3

    # Extract the necessary output

    # Calculate the work function

    # Add the workflow to the launchpad
    workflow = Workflow(fireworks=[relax_firework, dos_firework],
                        links_dict={relax_firework: [dos_firework]},
                        name=structure_file + " DOS calculation (Custodian)")

    LAUNCHPAD.add_wf(workflow)


VASP_BACKUP_FILES = {"INCAR", "KPOINTS", "POSCAR", "OUTCAR", "CONTCAR",
                     "OSZICAR", "vasprun.xml", "vasp.out", "std_err.txt"}

class QuotasErrorHandler(VaspErrorHandler):
    """
    Overwritten error handler for the issues we encounter often in slab
    calculations.

    """

    error_msgs = {
        "subspacematrix": ["WARNING: Sub-Space-Matrix is not hermitian in "
                           "DAV"],
        "edddav": ["Error EDDDAV: Call to ZHEGV failed"],
        "zbrent": ["ZBRENT: fatal internal in",
                   "ZBRENT: fatal error in bracketing"]
    }

    def __init__(self, output_filename="vasp.out", natoms_large_cell=100,
                 errors_subset_to_catch=None):
        """
        Initializes the handler with the output file to check.

        Args:
            output_filename (str): This is the file where the stdout for vasp
                is being redirected. The error messages that are checked are
                present in the stdout. Defaults to "vasp.out", which is the
                default redirect used by :class:`custodian.vasp.jobs.VaspJob`.
            natoms_large_cell (int): Number of atoms threshold to treat cell
                as large. Affects the correction of certain errors. Defaults to
                100.
            errors_subset_to_detect (list): A subset of errors to catch. The
                default is None, which means all supported errors are detected.
                Use this to only catch only a subset of supported errors.
                E.g., ["eddrrm", "zheev"] will only catch the eddrmm and zheev
                errors, and not others. If you wish to only excluded one or
                two of the errors, you can create this list by the following
                lines:

                ```
                subset = list(QuotasErrorHandler.error_msgs.keys())
                subset.pop("eddrrm")

                handler = QuotasErrorHandler(errors_subset_to_catch=subset)
                ```
        """
        super(QuotasErrorHandler, self).__init__(output_filename=output_filename,
                                      natoms_large_cell=natoms_large_cell)
        self.errors_subset_to_catch = errors_subset_to_catch or \
            list(QuotasErrorHandler.error_msgs.keys())

    def correct(self):
        backup(VASP_BACKUP_FILES | {self.output_filename})
        actions = []
        vi = VaspInput.from_directory(".")

        # If we encounter a DAV Sub-Space-Matrix error
        if "subspacematrix" or "edddav" in self.errors:
            # Switch to the CG algorithm
            actions.append(
                {"dict": "INCAR", "action": {"_set": {"ALGO": "All"}}}
            )

        # If we encounter a ZBRENT error
        if "zbrent" in self.errors:
            # Switch to Quasi-Newton algorithm
            actions.append({"dict": "INCAR",
                            "action": {"_set": {"IBRION": 1}}})
            # Move CONTCAR to POSCAR
            actions.append({"file": "CONTCAR",
                            "action": {"_file_copy": {"dest": "POSCAR"}}})

        VaspModder(vi=vi).apply_actions(actions)
        return {"errors": list(self.errors), "actions": actions}