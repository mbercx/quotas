# Encoding: UTF-8

import click

""" 
Command line interface for the quotas package.

"""

# This is only used to make '-h' a shorter way to access the CLI help
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

@click.group(context_settings=CONTEXT_SETTINGS)
def main():
    """
    A set of tools to help set up calculations and analyze results for the
    QUOTAS project.
    """
    pass

@main.group(context_settings=CONTEXT_SETTINGS)
def setup():
    """
    Set up all calculations for the input required in the QUOTAS model.
    """
    pass

@setup.command(context_settings=CONTEXT_SETTINGS)
@click.argument("miller_indices", nargs=1)
@click.argument("filename", nargs=1)
@click.option("--vacuum", "-V", default=float(15),
              help="Minimum thickness of the vacuum layer.")
@click.option("--thickness", "-t", default=20,
              help="Minimum thickness of the slab, in Angstroms.")
@click.option("--fix_part", "-f", default="center",
              help="Part of the slab to fix in the geometry optimization.")
@click.option("--fix_thickness", "-b", default=8,
              help="Number of layers fixed as bulk in the geometry "
                   "optimization.")
@click.option("-verbose", "-v", is_flag=True)
def slab(miller_indices, filename, vacuum, thickness, fix_part, fix_thickness,
         verbose):
    """
    Set up all the calculations for a specific surface of a structure.
    """
    from quotas.cli.commands.setup import slab_setup

    #TODO Add checks for the miller_indices
    miller_indices = [int(number) for number in miller_indices]

    slab_setup(filename=filename,
               miller_indices=miller_indices,
               thickness=thickness,
               vacuum=vacuum,
               fix_part=fix_part,
               fix_thickness=fix_thickness,
               verbose=verbose)