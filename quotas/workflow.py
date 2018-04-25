# Encoding: UTF-8

import os

from fireworks import FireTaskBase, Firework, LaunchPad, ScriptTask, \
    TemplateWriterTask, PyTask, \
    FileTransferTask, FWorker, Workflow, FWAction
from fireworks.core.rocket_launcher import launch_rocket, rapidfire
from fw_tutorials.firetask.addition_task import AdditionTask

"""
Workflow setup for the quotas calculations.

"""

TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            "templates")

class slabRelaxTask(FireTaskBase):
    """
    FireTask that sets up a slab geometry optimization.

    """
    _fw_name = "Slab Relaxation"

    def run_task(self, fw_spec):
        structure_file = fw_spec["structure_file"]
        fix_part = fw_spec["fix_part"]
        fix_thickness = fw_spec["fix_thickness"]
        is_metal = fw_spec["is_metal"]

        from quotas.cli.commands.slab import relax

        relax(structure_file, fix_part, fix_thickness, is_metal)

        return FWAction()


def dos_workflow(structure_file, fix_part, fix_thickness, is_metal):
    """
    Finally time for the real deal. I want to calculate the DOS and
    workfunction based on a structure file.


    Returns:

    """
    ## Various steps lined up. These should all be FireTasks within FireWorks


    # Set up the Launchpad for the workflow
    launchpad = LaunchPad(host="ds135179.mlab.com", port=35179, name="quotas",
                          username="mbercx", password="quotastests")

    ## FireWork 1

    # Set up the geometry optimization from the structure file
    setup_relax = PyTask(func="quotas.cli.commands.slab.relax",
                         kwargs=[structure_file,
                                 fix_part,
                                 fix_thickness,
                                 is_metal])

    # Set up the job script
    # TODO Allow scripts for various clusters

    job_script = TemplateWriterTask(
        {"template_file":os.path.join(TEMPLATE_DIR, "job_leibniz.sh"),
         "context":{"name":structure_file[:6] + "_rel", "nodes": 4},
         "output_file":"job_leibniz.sh"}
    )

    # Run the jobscript.
    job_submission = ScriptTask.from_str("msub job_leibniz.sh")

    fw = Firework([setup_relax, job_script, job_submission])

    # Launch the workflow
    launchpad.add_wf(fw)
    launch_rocket(launchpad)

    # -----> Here we would add a check to see if the job completed
    # successfully. If not, we can add another FireWork that makes the
    # necessary adjustments and restarts the calculation.

    ## Firework 2

    # Extract the necessary output from the geometry optimization, such as
    # geometry (duh) and the magnetic moments, MORE?
    # Create a custom FireTask from the corresponding quotas cli method?

    # Set up the job script

    # Run the job script

    # ----> Here we would add a check...

    ## Firework 3

    # Extract the necessary output

    # Calculate the work function


# Some tutorial-based tests:

def launch_test():
    """
    Small script that simply tests launching rockets using fireworks.
    """

    # Set up the launchpad, i.e. the connection with the mongo DB
    launchpad = LaunchPad(host="ds135179.mlab.com", port=35179, name="quotas",
                          username="mbercx", password="quotastests")
    launchpad.reset("", require_password=False)

    # Set up the firetask/work
    firetask = ScriptTask.from_str("echo 'Howdy! Your test was successful!'")
    firework = Firework(firetask)

    # Store workflow and launch
    launchpad.add_wf(firework)
    launch_rocket(launchpad)


def firework_test():
    """
    Script that tests the use of a Firework, consisting of various Firetasks.
    """

    # Set up the launchpad and reset it
    launchpad = LaunchPad(host="ds135179.mlab.com", port=35179, name="quotas",
                          username="mbercx", password="quotastests")
    launchpad.reset("", require_password=False)

    # Create the various firetasks for the firework
    firetask1 = TemplateWriterTask(
        {"context": {"opt1": 5.0, "opt2": "fast method"}, "template_file":
            "simple_template.txt", "output_file": "inputs.txt"}
    )
    firetask2 = ScriptTask.from_str("wc -w < inputs.txt > words.txt")
    firetask3 = FileTransferTask({"files": [{"src": "words.txt", "dest":
        "~/words.txt"}], "mode": "copy"})

    fw = Firework([firetask1, firetask2, firetask3])

    # Store workflow on the mongoDB and launch it
    launchpad.add_wf(fw)
    launch_rocket(launchpad)


def firetask_test():
    """
    It's a script that tests our very own firetask, dummy.

    """

    # Set up the launchpad and reset it
    launchpad = LaunchPad(host="ds135179.mlab.com", port=35179, name="quotas",
                          username="mbercx", password="quotastests")
    launchpad.reset("", require_password=False)

    # create the Firework consisting of a custom "Addition" task
    firework = Firework(AdditionTask(), spec={"input_array": [1, 2]})

    # store workflow and launch it locally
    launchpad.add_wf(firework)
    launch_rocket(launchpad, FWorker())


def workflow_test():
    """
    Take a guess.

    """

    # Set up the launchpad and reset it
    launchpad = LaunchPad(host="ds135179.mlab.com", port=35179, name="quotas",
                          username="mbercx", password="quotastests")
    launchpad.reset("", require_password=False)

    # define four individual FireWorks used in the Workflow
    task1 = ScriptTask.from_str('echo "Ingrid is the CEO."')
    task2 = ScriptTask.from_str('echo "Jill is a manager."')
    task3 = ScriptTask.from_str('echo "Jack is a manager."')
    task4 = ScriptTask.from_str('echo "Kip is an intern."')

    fw1 = Firework(task1)
    fw2 = Firework(task2)
    fw3 = Firework(task3)
    fw4 = Firework(task4)

    # assemble Workflow from FireWorks and their connections by id
    workflow = Workflow([fw1, fw2, fw3, fw4],
                        {fw1: [fw2, fw3], fw2: [fw4], fw3: [fw4]})

    # First give the list of fireworks, then a dictionary of 'family tree'
    # or links.

    # store workflow and launch it locally
    launchpad.add_wf(workflow)
