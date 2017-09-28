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

@main.command(context_settings=CONTEXT_SETTINGS)
@click.argument("filename")
def setup(filename):
    """
    Set up all the calculations for the input required for the QUOTAS model,
    based on the bulk structure.
    """
    from quotas.cli.commands.setup import slab_setup

    slab_setup(filename=filename)