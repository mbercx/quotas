# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import click

""" 
Command line interface for the quotas package.

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.2"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Apr 2018"

# This is used to make '-h' a shorter option to access the CLI help
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
@click.option("--vacuum", "-V", default=float(20),
              help="Minimum thickness of the vacuum layer, in angstrom. "
                   "Defaults to 20 angstrom.")
@click.option("--thickness", "-t", default=20,
              help="Minimum thickness of the slab, in angstrom. Defaults to "
                   "20 angstrom.")
@click.option("--write_cif", "-w", is_flag=True)
@click.option("--verbose", "-v", is_flag=True)
def setup(bulk_file, miller_indices, vacuum, thickness, write_cif, verbose):
    """
    Set up slabs from the bulk geometry.
    """
    from quotas.cli.commands.slab import setup

    #TODO Add checks for the miller_indices
    miller_indices = [int(number) for number in miller_indices]

    setup(bulk_file, miller_indices, thickness, vacuum, write_cif, verbose)


@slab.command(context_settings=CONTEXT_SETTINGS)
@click.argument("slab_file", nargs=1)
@click.option("--fix_part", "-f", default="center",
              help="Part of the slab to fix in the geometry optimization. "
                   "Defaults to 'center', which is currently the only "
                   "option. So it's not much of an option, really.")
@click.option("--fix_thickness", "-b", default=8,
              help="Number of layers fixed as bulk in the geometry "
                   "optimization. Defaults to 8.")
@click.option("--is_metal", "-m", is_flag=True)
@click.option("--verbose", "-v", is_flag=True)
def relax(slab_file, fix_part, fix_thickness, is_metal, verbose):
    """
    Set up the geometry optimization.
    """
    from quotas.cli.commands.slab import relax

    relax(slab_file, fix_part, fix_thickness, is_metal, verbose)


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
@click.option("--k_product", "-k", default=80,
              help="Product of the number of k-point along ")
@click.option("--hse", "-H", is_flag=True)
def dos(relax_dir, k_product, hse):
    """
    Set up the Density of states calculation.
    """
    from quotas.cli.commands.slab import dos

    dos(relax_dir=relax_dir,
        k_product=k_product,
        hse_calc=hse)


@main.group(context_settings=CONTEXT_SETTINGS)
def bulk():
    """
    Set up bulk calculations.
    """
    pass


@bulk.command(context_settings=CONTEXT_SETTINGS)
@click.argument("bulk_file", nargs=1)
@click.option("--is_metal", "-m", is_flag=True)
@click.option("--hse_calc", "-H", is_flag=True)
@click.option("--verbose", "-v", is_flag=True)
def relax(bulk_file, is_metal, hse_calc, verbose):
    """
    Set up the geometry optimization.
    """
    from quotas.cli.commands.bulk import relax

    relax(structure_file=bulk_file,
          is_metal=is_metal,
          hse_calc=hse_calc,
          verbose=verbose)


@bulk.command(context_settings=CONTEXT_SETTINGS)
@click.argument("relax_dir", nargs=1)
@click.option("--k_product", "-k", default=80,
              help="HELP")
@click.option("--hse_calc", "-H", is_flag=True)
def dos(relax_dir, k_product, hse_calc):
    """
    Set up the Density of states calculation.
    """
    from quotas.cli.commands.bulk import dos

    dos(relax_dir=relax_dir,
        k_product=k_product,
        hse_calc=hse_calc)


@bulk.command(context_settings=CONTEXT_SETTINGS)
@click.argument("relax_dir", nargs=1)
@click.option("--k_product", "-k", default=80)
@click.option("--hse_calc", "-H", is_flag=True)
@click.option("--is_metal", "-m", is_flag=True)
@click.option("--verbose", "-v", is_flag=True)
def optics(relax_dir, k_product, hse_calc, is_metal, verbose):
    """
    Set up a dielectric function calculation.
    """
    from quotas.cli.commands.bulk import optics

    optics(relax_dir=relax_dir,
           k_product=k_product,
           hse_calc=hse_calc,
           is_metal=is_metal,
           verbose=verbose)


@main.group(context_settings=CONTEXT_SETTINGS)
def analyze():
    """
    Analysis commands.
    """
    pass


@analyze.command(context_settings=CONTEXT_SETTINGS)
@click.argument("directory", nargs=1)
@click.option("--plot", "-p", is_flag=True)
def wf(directory, plot):
    """
    Calculate the work function from the local potential.
    """
    from quotas.cli.commands.analyze import wf

    wf(directory=directory,
       plot_potential=plot)


@analyze.command(context_settings=CONTEXT_SETTINGS)
@click.argument("directory", nargs=1)
def bandgap(directory):
    """
    Calculate the band gap from the electron energy levels.
    """
    from quotas.cli.commands.analyze import bandgap

    bandgap(directory=directory)

@analyze.command(context_settings=CONTEXT_SETTINGS)
@click.argument("vasprun", nargs=1)
def pdos(vasprun):
    """
    Calculate the band gap from the electron energy levels.
    """
    from quotas.cli.commands.analyze import pdos

    pdos(vasprun=vasprun)


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


@util.command(context_settings=CONTEXT_SETTINGS)
@click.option("--directory", "-d", default=".")
def process(directory):
    """
    Process the data into a .json file.
    """
    from quotas.cli.commands.util import process_output

    process_output(directory=directory)


@util.command(context_settings=CONTEXT_SETTINGS)
@click.argument("alloy_element", nargs=1)
@click.argument("structure_file", nargs=1)
@click.option("--mixing_ratio", "-M", default=0.5)
@click.option("--number_random", "-N", default=10)
def alloy(alloy_element, structure_file, mixing_ratio, number_random):
    """
    Set up a bunch of random structures to study allows.
    """
    from quotas.cli.commands.util import generate_alloy

    generate_alloy(root_structure_file=structure_file,
                   alloy_element=alloy_element,
                   mixing_ratio=mixing_ratio,
                   number_random=number_random)


@main.group(context_settings=CONTEXT_SETTINGS)
def workflow():
    """
    Set up workflows for quotas calculations.
    """
    pass


@workflow.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--fix_part", "-f", default="center",
              help="Part of the slab to fix in the geometry optimization. "
                   "Defaults to 'center', which is currently the only "
                   "option. So it's not much of an option, really.")
@click.option("--fix_thickness", "-b", default=8,
              help="Number of layers fixed as bulk in the geometry ")
@click.option("--is_metal", "-m", is_flag=True,
              help="Option that indicates that the material of the slab is "
                   "metallic. This is used to change the smearing settings "
                   "in the geometry optimization.")
@click.option("--k_product", "-k", default=80,
              help="Determines the density of the k-mesh in the density of "
                   "states calculation. k_product represents the product of "
                   "the number of k-points corresponding to a certain "
                   "lattice vector with the length of that lattice vector.")
@click.option("--in_custodian", "-c", is_flag=True,
              help="Specify that the VASP calculations of the workflow "
                   "should be run in a custodian.")
def dos(structure_file, fix_part, fix_thickness, is_metal, k_product,
        in_custodian):
    """
    Set up a DOS workflow script.
    """
    from quotas.workflow import dos_workflow

    dos_workflow(structure_file=structure_file,
                 fix_part=fix_part,
                 fix_thickness=fix_thickness,
                 is_metal=is_metal,
                 k_product=k_product,
                 in_custodian=in_custodian)


@workflow.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--is_metal", "-m", is_flag=True,
              help="Option that indicates that the material of the slab is "
                   "metallic. This is used to change the smearing settings "
                   "in the geometry optimization.")
@click.option("--hse_calc", "-H", is_flag=True,
              help="Sets the exchange-correlation functional of the "
                   "calculation to HSE06.")
@click.option("--k_product", "-k", default=80,
              help="Determines the density of the k-mesh in the density of "
                   "states calculation. k_product represents the product of "
                   "the number of k-points corresponding to a certain "
                   "lattice vector with the length of that lattice vector.")
@click.option("--in_custodian", "-c", is_flag=True,
              help="Specify that the VASP calculations of the workflow "
                   "should be run in a custodian.")
def optics(structure_file, is_metal, hse_calc, k_product, in_custodian):
    """
    Test the optics workflow script
    """
    from quotas.workflow import bulk_optics_workflow

    bulk_optics_workflow(structure_file=structure_file,
                         is_metal=is_metal,
                         hse_calc=hse_calc,
                         k_product=k_product,
                         in_custodian=in_custodian)


@main.group(context_settings=CONTEXT_SETTINGS)
def test():
    """
    Set of test scripts for various purposes.
    """
    pass


@test.command(context_settings=CONTEXT_SETTINGS)
@click.argument("directory", nargs=1)
def run_cust(directory):
    """
    Test the Custodian run.

    """
    from quotas.workflow import run_custodian

    run_custodian(directory)

@test.command(context_settings=CONTEXT_SETTINGS)
@click.argument("directory", nargs=1)
def check(directory):
    """
    Check a VASP run to see if it has completed succesfully.

    """
    from quotas.cli.commands.util import check_run

    check_run(directory=directory)


@test.command(context_settings=CONTEXT_SETTINGS)
@click.argument("structure_file", nargs=1)
@click.option("--fix_part", "-f", default="center",
              help="Part of the slab to fix in the geometry optimization. "
                   "Defaults to 'center', which is currently the only "
                   "option. So it's not much of an option, really.")
@click.option("--fix_thickness", "-b", default=8,
              help="Number of layers fixed as bulk in the geometry ")
@click.option("--is_metal", "-m", is_flag=True,
              help="Option that indicates that the material of the slab is "
                   "metallic. This is used to change the smearing settings "
                   "in the geometry optimization.")
@click.option("--k_product", "-k", default=80,
              help="Determines the density of the k-mesh in the density of "
                   "states calculation. k_product represents the product of "
                   "the number of k-points corresponding to a certain "
                   "lattice vector with the length of that lattice vector.")
def custodian(structure_file, fix_part, fix_thickness, is_metal, k_product):
    """
    Test the DOS workflow script
    """
    from quotas.workflow import test_custodian

    test_custodian(structure_file, fix_part, fix_thickness, is_metal, k_product)