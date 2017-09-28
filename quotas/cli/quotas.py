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
def setup(filename):
    """
    Set up all calculations for the input required in the QUOTAS model.
    """
    pass

@setup.command(context_settings=CONTEXT_SETTINGS)
@click.argument("miller_indices")
@click.argument("filename")
@click.option("--vacuum", "-v", default=15,
              help="Minimum thickness of the vacuum layer.")
@click.option("--thickness", "-t", default=20,
              help="Minimum thickness of the slab, in Angstroms.")
@click.option("--fix", "-f", default="center",
              help="Part of the slab to fix in the geometry optimization.")
def slab(miller_indices, filename):
    """
    Set up all the calculations for a specific surface of a structure.
    """
    from quotas.cli.commands.setup import slab_setup

    slab_setup(filename=filename,
               miller_indices=miller_indices)