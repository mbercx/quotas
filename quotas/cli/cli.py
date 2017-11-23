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


@main.group()
def slab():
    """
    Set up slab calculations.
    """
    pass


@slab.command(context_settings=CONTEXT_SETTINGS)
@click.argument("miller_indices", nargs=1)
@click.argument("bulk_file", nargs=1)
@click.option("--vacuum", "-V", default=float(15),
              help="Minimum thickness of the vacuum layer.")
@click.option("--thickness", "-t", default=20,
              help="Minimum thickness of the slab, in Angstroms.")
@click.option("--verbose", "-v", is_flag=True)
@click.option("--write_cif", "-c", is_flag=True)
def setup(bulk_file, miller_indices, vacuum, thickness, verbose, write_cif):
    """
    Set up slabs from the bulk geometry.
    """
    from quotas.cli.commands.slab import setup

    #TODO Add checks for the miller_indices
    miller_indices = [int(number) for number in miller_indices]

    setup(bulk_file, miller_indices, thickness, vacuum, verbose, write_cif)


@slab.command(context_settings=CONTEXT_SETTINGS)
@click.argument("slab_file", nargs=1)
@click.option("--fix_part", "-f", default="center",
              help="Part of the slab to fix in the geometry optimization.")
@click.option("--fix_thickness", "-b", default=8,
              help="Number of layers fixed as bulk in the geometry "
                   "optimization.")
@click.option("--verbose", "-v", is_flag=True)
def relax(slab_file, fix_part, fix_thickness, verbose):
    """
    Set up the geometry optimization.
    """
    from quotas.cli.commands.slab import relax

    relax(slab_file, fix_part, fix_thickness, verbose)


@slab.command(context_settings=CONTEXT_SETTINGS)
@click.argument("relax_dir", nargs=1)
@click.option("--k_product", "-k", default=50)
def wf(relax_dir, k_product):
    """
    Set up the work function calculation.
    """
    from quotas.cli.commands.slab import wf

    wf(relax_dir=relax_dir,
       k_product=k_product)


@slab.command(context_settings=CONTEXT_SETTINGS)
@click.argument("relax_dir", nargs=1)
@click.option("--k_product", "-k", default=80)
def dos(relax_dir, k_product):
    """
    Set up the Density of states calculation.
    """
    from quotas.cli.commands.slab import dos

    dos(relax_dir=relax_dir,
             k_product=k_product)


@main.group(context_settings=CONTEXT_SETTINGS)
def bulk():
    """
    Set up bulk calculations.
    """
    pass


@bulk.command(context_settings=CONTEXT_SETTINGS)
@click.argument("bulk_file", nargs=1)
@click.option("--verbose", "-v", is_flag=True)
def relax(bulk_file, verbose):
    """
    Set up the geometry optimization.
    """
    from quotas.cli.commands.bulk import relax

    relax(bulk_file=bulk_file,
          verbose=verbose)


@bulk.command(context_settings=CONTEXT_SETTINGS)
@click.argument("relax_dir", nargs=1)
@click.option("--k_product", "-k", default=80)
@click.option("--hse", "-H", is_flag=True)
def dos(relax_dir, k_product, hse):
    """
    Set up the Density of states calculation.
    """
    from quotas.cli.commands.bulk import dos

    dos(relax_dir=relax_dir,
        k_product=k_product,
        hse_calc=hse)


@bulk.command(context_settings=CONTEXT_SETTINGS)
def diel():
    """
    Set up a dielectric function calculation.
    """
    pass


@main.group(context_settings=CONTEXT_SETTINGS)
def analyze():
    """
    Analysis commands.
    """
    pass


@analyze.command(context_settings=CONTEXT_SETTINGS)
@click.argument("directory", nargs=1)
def wf(directory):
    """
    Calculate the work function from the local potential.
    """
    pass


@main.group(context_settings=CONTEXT_SETTINGS)
def util():
    """
    Utility commands.
    """
    pass


@util.command(context_settings=CONTEXT_SETTINGS)
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
    from quotas.cli.commands.util import kpar

    kpar(directory=directory,
         max_kpar=max_kpar,
         add_kpar=add_kpar)


@util.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".")
def nkp(directory):
    """
    Quickly find the number of kpoints based on the input files in the current
    directory.
    """
    from quotas.cli.commands.util import nkp

    nkp(directory=directory)