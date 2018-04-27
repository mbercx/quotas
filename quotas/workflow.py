# Encoding: UTF-8

import os
import subprocess

from fireworks import FireTaskBase, Firework, LaunchPad, PyTask, FWorker, \
    Workflow, FWAction
from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter

"""
Workflow setup for the quotas calculations.

"""

# Directory with templates for the Template
TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            "templates")
VASP_RUN_SCRIPT = "/user/antwerpen/202/vsc20248/local/scripts/job_workflow.sh"

# Set up the Launchpad for the workflows
LAUNCHPAD = LaunchPad(host="ds135179.mlab.com", port=35179, name="quotas",
                      username="mbercx", password="quotastests")

# Set up the FireWorker
FIREWORKER = FWorker(name="leibniz")

# Set up the queue adapter
QUEUE_ADAPTER = CommonAdapter.from_dict(
    {"_fw_q_type": "PBS",
     "rocket_launch": "source ~/local/envs/pymatgen.env; "
                      "rlaunch singleshot",
     "nnodes": "1",
     "ppnode": "28",
     "walltime": "72:00:00",
     "queue": "batch",
     "job_name": "test",
     "logdir": "/user/antwerpen/202/vsc20248", }
)


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


def dos_workflow(structure_file, fix_part, fix_thickness, is_metal, k_product):
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
    run_relax = PyTask(func="quotas.workflow.run_vasp",
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
    run_dos = PyTask(func="quotas.workflow.run_vasp",
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
                        name=structure_file + " DOS calculation")

    LAUNCHPAD.add_wf(workflow)


def bulk_optics_workflow(structure_file, is_metal, hse_calc, k_product):
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

    # Run VASP
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

    # Run VASP
    run_optics = PyTask(func="quotas.workflow.run_vasp",
                        inputs=["optics_dir"])

    optics_firework = Firework(tasks=[setup_optics, run_optics],
                               name="Optics calculation")

    # Add the workflow to the launchpad
    workflow = Workflow(fireworks=[relax_firework, optics_firework],
                        links_dict={relax_firework: [optics_firework]},
                        name=structure_file + " Optics calculation")

    LAUNCHPAD.add_wf(workflow)
