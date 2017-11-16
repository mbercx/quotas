# Encoding: UTF-8

import click

""" 
Command line interface for the quotas package.

"""

# This is used to make '-h' a shorter way to access the CLI help
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

@click.group(context_settings=CONTEXT_SETTINGS)
def main():
    """
    A set of tools to help set up calculations and analyze results for the
    QUOTAS project.
    """
    pass

@main.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".")
@click.option("--max_kpar", "-m", default=30,
              help="Maximum value for KPAR, usually determined by the "
                   "number of nodes at the user's disposal.")
@click.option("--add_kpar", "-a", is_flag=True,
              help="Add the KPAR tag to the INCAR file.")
def kpar(directory, max_kpar, add_kpar):
    """
    Find a suitable value for KPAR.
    """
    from quotas.cli.commands.setup import kpar

    kpar(directory=directory,
         max_kpar=max_kpar,
         add_kpar=add_kpar)


@main.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".")
def nkp(directory):
    """
    Quickly find the number of kpoints based on the input files in the current
    directory.
    """
    from quotas.cli.commands.setup import nkp

    nkp(directory=directory)


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
@click.option("--verbose", "-v", is_flag=True)
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

@setup.command(context_settings=CONTEXT_SETTINGS)
@click.argument("relax_dir", nargs=1)
@click.option("--k_product", "-k", default=50)
def wf(relax_dir, k_product):
    """
    Set up the work function calculation, based on the output of the geometry
    optimization.

    """
    from quotas.cli.commands.setup import work_function_calc

    work_function_calc(relax_dir=relax_dir,
                       k_product=k_product)


@setup.command(context_settings=CONTEXT_SETTINGS)
@click.argument("relax_dir", nargs=1)
@click.option("--k_product", "-k", default=80)
def dos(relax_dir, k_product):
    """
    Set up the DOS calculation, based on the output of the geometry
    optimization.

    """
    from quotas.cli.commands.setup import DOS_calc

    DOS_calc(relax_dir=relax_dir,
             k_product=k_product)