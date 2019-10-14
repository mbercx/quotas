# Encoding: UTF-8
# Copyright (c) Marnik Bercx, University of Antwerp
# Distributed under the terms of the MIT License

import numpy as np
import os, cmath
import copy
import string
import sys
import json
import warnings
import matplotlib.pyplot as plt

from scipy import constants
from scipy.constants import hbar, e, m_e

from fnmatch import fnmatch
from pymatgen import Lattice, PeriodicSite, Structure, Composition
from pymatgen.core.surface import SlabGenerator, Slab
from pymatgen.electronic_structure.dos import Dos
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Outcar, Locpot, Vasprun
from pymatgen.util.plotting import pretty_plot
from monty.json import MSONable, MontyDecoder, MontyEncoder
from monty.io import zopen

"""
A set of methods to aid in the setup of slab calculations. 

"""

__author__ = "Marnik Bercx"
__copyright__ = "Copyright 2018, Marnik Bercx, University of Antwerp"
__version__ = "0.2"
__maintainer__ = "Marnik Bercx"
__email__ = "marnik.bercx@uantwerpen.be"
__date__ = "Apr 2018"


class QSlab(Slab):
    """
    A Quotas version of the pymatgen.core.surface.Slab object. All methods can be
    inherited, but we need to add some convenience methods of our own.

    """

    def __init__(self, lattice, species, coords, miller_index,
                 oriented_unit_cell, shift, scale_factor, reorient_lattice=True,
                 validate_proximity=False, to_unit_cell=False,
                 reconstruction=None, coords_are_cartesian=False,
                 site_properties=None, energy=None):
        super(QSlab, self).__init__(lattice, species, coords, miller_index,
                                    oriented_unit_cell, shift, scale_factor,
                                    reorient_lattice, validate_proximity,
                                    to_unit_cell,
                                    reconstruction, coords_are_cartesian,
                                    site_properties, energy)

    @classmethod
    def from_file(cls, filename, primitive=False, sort=False, merge_tol=0.0):
        """
        Load the QSlab from a file.

        Args:
            filename (str): Path to file that contains the details of the slab.
            primitive (bool):
            sort (bool):
            merge_tol (float):

        Returns:

        """
        if fnmatch(os.path.basename(filename).lower(), "*.json"):
            with zopen(filename, "r") as file:
                return cls.from_str(file.read())
        else:
            raise NotImplementedError("Only .json files are currently supported.")

    def as_dict(self):
        d = super(Slab, self).as_dict()
        d["@module"] = self.__class__.__module__
        d["@class"] = self.__class__.__name__
        d["oriented_unit_cell"] = self.oriented_unit_cell.as_dict()
        d["miller_index"] = self.miller_index
        d["shift"] = self.shift
        d["scale_factor"] = MontyEncoder().default(self.scale_factor)
        d["reconstruction"] = self.reconstruction
        d["energy"] = self.energy
        return d

    @classmethod
    def from_dict(cls, d):
        lattice = Lattice.from_dict(d["lattice"])
        sites = [PeriodicSite.from_dict(sd, lattice) for sd in d["sites"]]
        s = Structure.from_sites(sites)

        return cls(
            lattice=lattice,
            species=s.species_and_occu, coords=s.frac_coords,
            miller_index=d["miller_index"],
            oriented_unit_cell=Structure.from_dict(d["oriented_unit_cell"]),
            shift=d["shift"],
            scale_factor=MontyDecoder().process_decoded(d["scale_factor"]),
            site_properties=s.site_properties, energy=d["energy"]
        )

    @classmethod
    def from_str(cls, input_string, fmt="json", primitive=False, sort=False,
                 merge_tol=0.0):

        if fmt == "json":
            return QSlab.from_slab(json.loads(input_string, cls=MontyDecoder))
        else:
            raise NotImplementedError("Currently only the json format is supported.")

    def to(self, fmt=None, filename=None, **kwargs):
        """
        Write the QSlab to a file.

        Args:
            fmt:
            filename:

        Returns:

        """
        if fmt == "json":
            with open(filename, "w") as file:
                file.write(self.to_json())
        else:
            super(QSlab, self).to(fmt, filename, **kwargs)

    def get_sorted_structure(self, key=None, reverse=False):
        """
        Get a sorted copy of the structure. The parameters have the same
        meaning as in list.sort. By default, sites are sorted by the
        electronegativity of the species. Note that Slab has to override this
        because of the different __init__ args.

        Args:
            key: Specifies a function of one argument that is used to extract
                a comparison key from each list element: key=str.lower. The
                default value is None (compare the elements directly).
            reverse (bool): If set to True, then the list elements are sorted
                as if each comparison were reversed.
        """
        sites = sorted(self, key=key, reverse=reverse)
        s = Structure.from_sites(sites)
        return self.__class__(s.lattice, s.species_and_occu, s.frac_coords,
                              self.miller_index, self.oriented_unit_cell, self.shift,
                              self.scale_factor, site_properties=s.site_properties,
                              reorient_lattice=self.reorient_lattice)

    def copy(self, site_properties=None, sanitize=False):
        """
        Convenience method to get a copy of the structure, with options to add
        site properties.

        Args:
            site_properties (dict): Properties to add or override. The
                properties are specified in the same way as the constructor,
                i.e., as a dict of the form {property: [values]}. The
                properties should be in the order of the *original* structure
                if you are performing sanitization.
            sanitize (bool): If True, this method will return a sanitized
                structure. Sanitization performs a few things: (i) The sites are
                sorted by electronegativity, (ii) a LLL lattice reduction is
                carried out to obtain a relatively orthogonalized cell,
                (iii) all fractional coords for sites are mapped into the
                unit cell.

        Returns:
            A copy of the Structure, with optionally new site_properties and
            optionally sanitized.
        """
        props = self.site_properties
        if site_properties:
            props.update(site_properties)
        return self.__class__(self.lattice, self.species_and_occu, self.frac_coords,
                              self.miller_index, self.oriented_unit_cell, self.shift,
                              self.scale_factor, site_properties=props,
                              reorient_lattice=self.reorient_lattice)

    @classmethod
    def from_slab(cls, slab):
        return cls(lattice=slab.lattice, species=slab.species,
                   coords=slab.frac_coords,
                   miller_index=slab.miller_index,
                   oriented_unit_cell=slab.oriented_unit_cell,
                   shift=slab.shift,
                   scale_factor=slab.scale_factor,
                   reconstruction=slab.reconstruction,
                   coords_are_cartesian=False,
                   site_properties=slab.site_properties, energy=slab.energy)

    def find_atomic_layers(self, layer_tol=2e-2):
        """
            Determines the atomic layers in the c-direction of the slab. Note that as
            long as a site is "close enough" to ONE other site of a layer (determined
            by the 'layer_tol' variable), it will be added to that layer.

            Another option would be to demand that the distance is smaller than
            'layer_tol' for ALL sites of the layer, but then the division in layers
            could depend on the order of the sites.

            Args:
                layer_tol (float): Tolerance for the maximum distance (in
                    angstrom) between layers for them to still correspond to the same
                    layer.

            Returns:
                (list) List of the atomic layers, sorted by their position in the c-direction.
                Each atomic layer is also represented by a list of sites.
            """

        atomic_layers = []

        # Get a unit vector perpendicular to the layers (i.e. the a and b lattice vector)
        m = self.lattice.matrix
        u = np.cross(m[0, :], m[1, :])
        u /= np.linalg.norm(u)
        c_proj = np.dot(u, m[2, :])

        for site in self.sites:

            is_in_layer = False

            # Check to see if the site is in a layer that is already in our list
            for layer in atomic_layers:

                # Compare the third fractional coordinate of the site with that of
                # the atoms in the considered layer
                for atom_site in layer.copy():

                    distance = abs(atom_site.frac_coords[2] - site.frac_coords[2]) * \
                               c_proj

                    if distance < layer_tol or abs(distance - c_proj) < layer_tol:
                        is_in_layer = True
                        layer.append(site)
                        break  # Break out of the loop, else the site is added
                        # multiple times

            # If the site is not found in any of the atomic layers, create a new
            # atomic layer
            if not is_in_layer:
                atomic_layers.append([site, ])

        # Sort the atomic layers
        atomic_layers.sort(key=lambda layer: layer[0].frac_coords[2])

        return atomic_layers

    def update_sites(self, directory, ignore_magmom=False):
        """
        Based on the CONTCAR and OUTCAR of a VASP calculation, update the
        site coordinates and magnetic moments of the slab.

        Args:
            directory (str): Directory in which the calculation output files (i.e.
                CONTCAR and OUTCAR) are stored.
            ignore_magmom (bool): Flag that indicates that the final magnetic
                moments of the calculation should be ignored.

        """
        new_slab = Structure.from_file(os.path.join(directory, "CONTCAR"))

        if ignore_magmom and "magmom" in self.site_properties.keys():
            new_slab.add_site_property("magmom", self.site_properties["magmom"])
        else:
            out = Outcar(os.path.join(directory, "OUTCAR"))
            if len(out.magnetization) == 0:
                if "magmom" in self.site_properties.keys():
                    warnings.warn("Outcar does not contain any magnetic moments! "
                                  "Keeping magnetic moments from initial slab.")
                    new_slab.add_site_property("magmom",
                                               self.site_properties["magmom"])
            else:
                new_slab.add_site_property("magmom",
                                           [site["tot"] for site in
                                            out.magnetization])

        # Update the lattice
        self._lattice = new_slab.lattice

        # Update the coordinates of the occupied sites.
        for i, site in enumerate(self):
            new_site = new_slab.sites[i]

            # Update the site coordinates
            self.replace(i, species=new_site.species,
                         coords=new_site.frac_coords,
                         properties=new_site.properties)


class QuotasCalculator(MSONable):
    """
    Calculator class for the QUOTAS project.

    """

    def __init__(self, total_dos, workfunction_data, dieltensor=None,
                 energy_spacing=1e-2, plasmon_parameters=None):
        """
        Initialize a QuotasCalculator.

        """

        self.total_dos = total_dos
        self.workfunction_data = workfunction_data
        self.dieltensor = dieltensor
        self.energy_spacing = energy_spacing
        self.plasmon_parameters = plasmon_parameters

        self.energies = np.arange(
            total_dos.energies.min(), total_dos.energies.max(), energy_spacing
        )
        self.dos = np.interp(self.energies, total_dos.energies,
                             sum(list(total_dos.densities.values())))

        self.escape_function = self.set_up_escape_function()

        self._bulk_plas_prob = None
        self._surf_plas_prob = None

        if dieltensor is not None:
            plasmon_parameters = plasmon_parameters or {"bulk": 0.095,
                                                        "surface": 1.6}
            self.set_up_plasmon_probabilities(
                bulk_parameter=plasmon_parameters["bulk"],
                surface_parameter=plasmon_parameters["surface"]
            )

    def set_up_escape_function(self):
        """
        Set up the escape function for the secondary electrons.

        TODO: This function should be optimized
        Right now this function has been largely copied in semantics from its MATLAB
        counterpart. I believe this can be more efficient.

        Returns:
            (numpy.ndarray): N-sized vector that represents the probability for
                an electron with a certain energy to escape.

        """

        conduction_energy = self.total_dos.get_cbm_vbm()[1]
        barrier = self.workfunction_data.vacuum_locpot - conduction_energy

        theta = np.arange(0, np.pi / 2, 0.01)
        return [np.trapz(np.array(
            [self.step_escape_probability(t, energy - conduction_energy, barrier)
             for t in theta]) * np.sin(theta), theta) / 2
                if energy > self.workfunction_data.vacuum_locpot else 0
                for energy in self.energies]

    def set_up_plasmon_probabilities(self, bulk_parameter, surface_parameter):
        """

        Args:
            surface_parameter:

        Returns:

        """
        self.plasmon_parameters = {"bulk": bulk_parameter,
                                   "surface": surface_parameter}
        bulk_loss_function = self.dieltensor.get_loss_function()
        conduction_energy = self.total_dos.get_cbm_vbm()[1]

        bulk_plas_prob = []
        diffs = np.diff(self.dieltensor.energies)
        max_index = 0
        total_loss_rate = 0

        for energy in self.energies:

            while self.dieltensor.energies[max_index + 1] \
                    < energy - conduction_energy:
                max_index += 1
                total_loss_rate += (bulk_loss_function[max_index] +
                                    bulk_loss_function[max_index - 1]) / 2 * \
                                   diffs[max_index - 1]

            if energy < conduction_energy or total_loss_rate == 0:
                bulk_plas_prob.append(np.zeros(bulk_loss_function.shape))

            else:
                loss_probabilities = bulk_loss_function / total_loss_rate * \
                                     (1 - np.exp(-total_loss_rate * bulk_parameter))
                loss_probabilities[self.dieltensor.energies
                                   > energy - conduction_energy] = 0
                bulk_plas_prob.append(loss_probabilities)

        bulk_plas_prob = np.array(bulk_plas_prob).swapaxes(0, 1)

        surface_loss_function = self.dieltensor.get_loss_function(surface=True)
        surface_plasmon_prob = surface_parameter * surface_loss_function / \
                               (surface_parameter * surface_loss_function + 1)

        surface_plasmon_prob = np.interp(self.energies, self.dieltensor.energies,
                                         surface_plasmon_prob)

        self._bulk_plas_prob = bulk_plas_prob
        self._surf_plas_prob = surface_plasmon_prob

    @property
    def bulk_plas_prob(self):
        return self._bulk_plas_prob

    @property
    def surf_plas_prob(self):
        return self._surf_plas_prob

    @property
    def occupied_states(self):
        occupied_states = self.dos.copy()
        occupied_states[self.energies > self.total_dos.efermi] = 0
        return occupied_states

    @property
    def empty_states(self):
        empty_states = self.dos.copy()
        empty_states[self.energies < self.total_dos.efermi] = 0
        return empty_states

    def calculate_yield(self, ion_energy, yield_convergence=1e-3):
        """

        Args:
            ion_energy:

        Returns:

        """

        excited_density = self.auger_neutralization(ion_energy)

        iteration_yield = 1
        yield_densities = []
        total_yields = []

        while iteration_yield > yield_convergence:

            if self.surf_plas_prob is not None:
                excited_density = self.bulk_plasmon_excitation(excited_density)

            yield_density, excited_density = self.electon_escape(excited_density)
            yield_densities.append(yield_density)

            excited_density = self.electron_scatter(excited_density)
            iteration_yield = np.trapz(yield_density, self.energies)
            total_yields.append(iteration_yield)

        vac_index = sum(self.energies < self.workfunction_data.vacuum_locpot)
        final_yield_density = sum(yield_densities)
        final_yield_density = final_yield_density[vac_index:]
        final_energies = self.energies[vac_index:] \
                         - self.workfunction_data.vacuum_locpot

        return final_energies, final_yield_density, sum(total_yields)

    def auger_neutralization(self, ion_energy):
        """
        Calculate the distribution of excited electrons after the Auger
        Neutralization of an incoming ion.

        Returns:

        """
        # Calculate the energy released upon neutralization
        neutralization_energy = np.roll(
            self.occupied_states,
            int((ion_energy - self.workfunction_data.vacuum_locpot) /
                self.energy_spacing)
        )

        if self.surf_plas_prob is not None:
            pre_norm = np.trapz(neutralization_energy, self.energies)
            neutralization_energy -= neutralization_energy * self.surf_plas_prob
            plasmon_fraction = np.trapz(neutralization_energy,
                                        self.energies) / pre_norm

        excited_density = np.convolve(
            self.occupied_states,
            neutralization_energy[sum(self.energies < 0):],
            "full"
        )

        excited_density = excited_density[:len(self.energies)]
        excited_density *= self.empty_states
        normalization = np.trapz(excited_density, self.energies)

        if self.surf_plas_prob is not None:
            normalization /= plasmon_fraction

        excited_density /= normalization

        return excited_density

    def electon_escape(self, excited_density):
        """
        Calculate the distribution of electrons that have escaped.

        Returns:

        """
        yield_density = excited_density * self.escape_function
        new_excited_density = excited_density - yield_density

        return yield_density, new_excited_density

    def electron_scatter(self, excited_density):
        """
        Calculate the distribution of excited electrons after scattering.

        Returns:

        """
        # Remove electrons which are below the vacuum level
        conduction_dos = self.empty_states
        excited_density[self.energies < self.workfunction_data.vacuum_locpot] = 0

        electrons_left = np.trapz(excited_density, self.energies)

        energy_distribution = np.convolve(self.occupied_states, excited_density,
                                          "full")

        scatter_transform = []
        for i in range(len(self.energies)):
            scatter_transform.append(energy_distribution[i:i + len(self.energies)])
        scatter_transform = np.array(scatter_transform)

        conduction_matrix = np.matmul(
            np.reshape(conduction_dos, (-1, 1)),
            np.reshape(conduction_dos, (-1, 1)).transpose()
        )

        scatter_transform *= conduction_matrix

        scatter_distribution = np.trapz(scatter_transform, self.energies, 0)
        scatter_distribution *= 2 * electrons_left / np.trapz(
            scatter_distribution, self.energies)

        return scatter_distribution

    def bulk_plasmon_excitation(self, excited_density):
        """
        Allow for a bulk plasmon excitation event for the excited density of
        electrons.

        Args:
            excited_density (np.ndarray):

        Returns:

        """
        plasmon_density = []
        for row in self.bulk_plas_prob:
            plasmon_density.append(row * excited_density)
        plasmon_density = np.array(plasmon_density)
        plasmon_loss = np.trapz(plasmon_density, self.dieltensor.energies, axis=0)

        plasmon_final = []
        zero_row = np.zeros(self.energies.shape)

        for plasmon_energy, row in zip(self.dieltensor.energies, plasmon_density):

            if np.all(row == 0):
                plasmon_final.append(zero_row)
            else:
                number_plasmons = np.trapz(row, self.energies)
                e_after_energy_loss = np.roll(row, -int(plasmon_energy /
                                                        self.energy_spacing))
                e_from_plasmon_decay = np.roll(
                    self.occupied_states, int(plasmon_energy / self.energy_spacing))
                e_from_plasmon_decay[self.empty_states < 0] = 0
                e_from_plasmon_decay *= number_plasmons / np.trapz(
                    e_from_plasmon_decay, self.energies)
                plasmon_final.append(
                    e_after_energy_loss + e_from_plasmon_decay
                )

        plasmon_final = np.array(plasmon_final)
        final_density = np.trapz(plasmon_final, self.dieltensor.energies, axis=0)

        return excited_density - plasmon_loss + final_density

    def as_dict(self):
        """


        Returns:

        """
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "total_dos": self.total_dos.as_dict(),
            "workfunction_data": self.workfunction_data.as_dict(),
            "dieltensor": self.dieltensor.as_dict(),
            "energy_spacing": self.energy_spacing,
            "plasmon_parameters": self.plasmon_parameters
        }
        return d

    def to(self, filename):
        """
        Write the QuotasCalculator to a JSON file.

        Args:
            filename (str): Path to the file.

        Returns:
            None

        """
        with zopen(filename, "w") as f:
            f.write(self.to_json())


    @classmethod
    def from_dict(cls, d):
        """

        Args:
            d:

        Returns:

        """
        return cls(
            total_dos=Dos.from_dict(d["total_dos"]),
            workfunction_data=WorkFunctionData.from_dict(d["workfunction_data"]),
            dieltensor=DielTensor.from_dict(d["dieltensor"]),
            energy_spacing=d["energy_spacing"],
            plasmon_parameters=d["plasmon_parameters"]
        )

    @classmethod
    def from_file(cls, filename, energy_spacing_override=None):
        """
        Initialize a QuotasCalculator instance from a file.

        Args:
            filename (str): Path to JSON file.
        Returns:

        """

        # JSON format
        with zopen(filename, "r") as f:
            d = json.loads(f.read())
            if energy_spacing_override is not None:
                d["energy_spacing"] = energy_spacing_override
            return cls.from_dict(d)

    @staticmethod
    def step_escape_probability(angle, energy, barrier):
        """
        Calculate the escape probability of a step function barrier.

        Args:
            angle:
            energy:
            barrier:

        Returns:

        """
        if energy < barrier:
            return 0

        k = (2 * energy * m_e / (hbar / e) ** 2) ** (1 / 2)

        k_perp = k * np.cos(angle)
        energy_perp = k_perp ** 2 * (hbar / e) ** 2 / 2 / m_e

        if energy_perp > barrier:
            p_perp = (2 * m_e * (energy_perp - barrier)) ** (1 / 2) / (hbar / e)
        else:
            p_perp = 0

        return 4 * k_perp * p_perp / (k_perp + p_perp) ** 2


# Stolen from pymatgen.analysis.surface_analysis
# Idea is to write a nice postprocessing tool extracting the necessary
# output from the work function and DOS calculations.

class WorkFunctionData(MSONable):
    """
    A class that contains the necessary data to calculate the work function
    of a slab, as well as plot the average local potential along the line
    perpendicular to the surface.

    Shamelessly copied from pymatgen.analysis.surface_analysis, in order to
    create a more lean class which I can store in a json file and push to
    the dataserver.

    """

    def __init__(self, poscar, locpot_along_c, efermi, shift=0):
        """
        Initializes the WorkFunctionData class.

        Args:
            poscar (pymatgen.io.vasp.outputs.Poscar): Poscar vasp output object
            locpot_along_c (numpy.ndarray): Average potential along the c axis.
            efermi (float): Fermi level, i.e. the energy of the highest
                occupied state.
            shift (float): Parameter to translate the slab (and
                therefore the vacuum) of the slab structure, thereby
                translating the plot along the x axis.
        """

        self.poscar = poscar

        # properties that can be shifted
        slab = poscar.structure.copy()
        slab.translate_sites([i for i, site in enumerate(slab)], [0, 0, shift])
        self.slab = slab
        self.sorted_sites = sorted(self.slab,
                                   key=lambda site: site.frac_coords[2])

        # Get the plot points between 0 and c
        # increments of the number of locpot points
        self.along_c = np.linspace(0, 1, num=len(locpot_along_c))
        locpot_along_c_mid, locpot_end, locpot_start = [], [], []
        for i, s in enumerate(self.along_c):
            j = s + shift
            if j > 1:
                locpot_start.append(locpot_along_c[i])
            elif j < 0:
                locpot_end.append(locpot_along_c[i])
            else:
                locpot_along_c_mid.append(locpot_along_c[i])
        self.locpot_along_c = locpot_start + locpot_along_c_mid + locpot_end

        # get the average of the signal in the bulk-like region of the
        # slab, i.e. the average of the oscillating region. This gives
        # a rough appr. of the potential in the interior of the slab
        bulk_p = [p for i, p in enumerate(self.locpot_along_c) if \
                  self.sorted_sites[-1].frac_coords[2] > self.along_c[i] \
                  > self.sorted_sites[0].frac_coords[2]]
        self.ave_bulk_p = np.mean(bulk_p)

        # shift independent quantities
        self.efermi = efermi
        self.vacuum_locpot = max(self.locpot_along_c)
        # get the work function
        self.work_function = self.vacuum_locpot - self.efermi
        # for setting ylim and annotating
        self.ave_locpot = (self.vacuum_locpot - min(self.locpot_along_c)) / 2

    def get_locpot_along_slab_plot(self, label_energies=True, plt=None,
                                   label_fontsize=10):
        """
        Returns a plot of the local potential (eV) vs the
            position along the c axis of the slab model (Ang)

        Args:
            label_energies (bool): Whether to label relevant energy
                quantities such as the work function, Fermi energy,
                vacuum locpot, bulk-like locpot
            plt (plt): Matplotlib pylab object
            label_fontsize (float): Fontsize of labels

        Returns plt of the locpot vs c axis
        """

        plt = pretty_plot() if not plt else plt
        ax = list(plt.subplots())[1]
        ax.spines['top'].set_visible(False)

        # plot the raw locpot signal along c
        plt.plot(self.along_c, self.locpot_along_c, 'b--')

        # Get the local averaged signal of the locpot along c
        xg, yg = [], []
        for i, p in enumerate(self.locpot_along_c):
            # average signal is just the bulk-like potential when in the slab region
            if p < self.ave_bulk_p \
                    or self.sorted_sites[-1].frac_coords[2] >= self.along_c[i] \
                    >= self.sorted_sites[0].frac_coords[2]:
                yg.append(self.ave_bulk_p)
                xg.append(self.along_c[i])
            else:
                yg.append(p)
                xg.append(self.along_c[i])
        xg, yg = zip(*sorted(zip(xg, yg)))
        plt.plot(xg, yg, 'r', linewidth=2.5, zorder=-1)

        # make it look nice
        if label_energies:
            plt = self.get_labels(plt, label_fontsize=label_fontsize)
        plt.xlim([0, 1])
        plt.ylim([min(self.locpot_along_c),
                  self.vacuum_locpot + self.ave_locpot * 0.2])
        plt.xlabel(r"Fractional coordinates ($\hat{c}$)", fontsize=25)
        plt.xticks(fontsize=15, rotation=45)
        plt.ylabel(r"Energy (eV)", fontsize=25)
        plt.yticks(fontsize=15)

        return plt

    def get_labels(self, plt, label_fontsize=10):
        """
        Handles the optional labelling of the plot with relevant quantities
        Args:
            plt (plt): Plot of the locpot vs c axis
            label_fontsize (float): Fontsize of labels
        Returns Labelled plt
        """

        maxc = self.sorted_sites[-1].frac_coords[2]
        minc = self.sorted_sites[0].frac_coords[2]
        # determine whether most of the vacuum is to
        # the left or right for labelling purposes
        vleft = [i for i in self.along_c if i <= minc]
        vright = [i for i in self.along_c if i >= maxc]
        if max(vleft) - min(vleft) > max(vright) - min(vright):
            label_in_vac = (max(vleft) - min(vleft)) / 2
        else:
            label_in_vac = (max(vright) - min(vright)) / 2 + maxc

        # label the vacuum locpot
        label_in_bulk = maxc - (maxc - minc) / 2
        plt.plot([0, 1], [self.vacuum_locpot] * 2, 'b--', zorder=-5,
                 linewidth=1)
        xy = [label_in_bulk, self.vacuum_locpot + self.ave_locpot * 0.05]
        plt.annotate(r"$V_{vac}=%.2f$" % (self.vacuum_locpot), xy=xy,
                     xytext=xy, color='b', fontsize=label_fontsize)

        # label the fermi energy
        plt.plot([0, 1], [self.efermi] * 2, 'g--',
                 zorder=-5, linewidth=3)
        xy = [label_in_bulk, self.efermi + self.ave_locpot * 0.05]
        plt.annotate(r"$E_F=%.2f$" % (self.efermi), xytext=xy,
                     xy=xy, fontsize=label_fontsize, color='g')

        # label the bulk-like locpot
        plt.plot([0, 1], [self.ave_bulk_p] * 2, 'r--', linewidth=1., zorder=-1)
        print(label_in_vac)
        xy = [label_in_vac, self.ave_bulk_p + self.ave_locpot * 0.05]
        plt.annotate(r"$V^{interior}_{slab}=%.2f$" % (self.ave_bulk_p),
                     xy=xy, xytext=xy, color='r', fontsize=label_fontsize)

        # label the work function as a barrier
        plt.plot([label_in_vac] * 2, [self.efermi, self.vacuum_locpot],
                 'k--', zorder=-5, linewidth=2)
        xy = [label_in_vac, self.efermi + self.ave_locpot * 0.05]
        plt.annotate(r"$\Phi=%.2f$" % (self.work_function),
                     xy=xy, xytext=xy, fontsize=label_fontsize)

        return plt

    def as_dict(self):
        """

        Returns:

        """

        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "poscar": self.poscar.as_dict(),
             "locpot_along_c": self.locpot_along_c,
             "efermi": self.efermi}

        return d

    @classmethod
    def from_dict(cls, d):
        """

        Args:
            d:

        Returns:

        """

        return cls(Poscar.from_dict(d["poscar"]), d["locpot_along_c"],
                   d["efermi"])

    @classmethod
    def from_output(cls, poscar_filename, locpot_filename, outcar_filename,
                    shift=0):
        poscar = Poscar.from_file(poscar_filename)
        locpot = Locpot.from_file(locpot_filename)
        outcar = Outcar(outcar_filename)

        wf = WorkFunctionAnalyzer(poscar, locpot, outcar, shift=shift)

        return cls(poscar, np.array(wf.locpot_along_c), wf.efermi,
                   shift=shift)

    @classmethod
    def from_dir(cls, directory):
        return cls.from_output(
            poscar_filename=os.path.join(directory, "CONTCAR"),
            locpot_filename=os.path.join(directory, "LOCPOT"),
            outcar_filename=os.path.join(directory, "OUTCAR")
        )

    @classmethod
    def from_file(cls, filename):
        """»

        Args:
            filename:

        Returns:

        """
        with open(filename, "r") as file:
            return json.load(file, cls=MontyDecoder)

    def to(self, filename):
        """

        Args:
            filename:

        Returns:

        """
        with open(filename, "w") as file:
            json.dump(self.as_dict(), file, cls=MontyEncoder)


class WorkFunctionAnalyzer(MSONable):
    """
    A class that post processes a task document of a vasp calculation (from
        using drone.assimilate). Can calculate work function from the vasp
        calculations and plot the potential along the c axis. This class
        assumes that LVTOT=True (i.e. the LOCPOT file was generated) for a
        slab calculation and it was insert into the task document along with
        the other outputs.

        Own version, to make it JSON serializable

    .. attribute:: efermi

        The Fermi energy

    .. attribute:: locpot_along_c

        Local potential in eV along points along the  axis

    .. attribute:: vacuum_locpot

        The maximum local potential along the c direction for
            the slab model, ie the potential at the vacuum

    .. attribute:: work_function

        The minimum energy needed to move an electron from the
            surface to infinity. Defined as the difference between
            the potential at the vacuum and the Fermi energy.

    .. attribute:: slab

        The slab structure model

    .. attribute:: along_c

        Points along the c direction with same
            increments as the locpot in the c axis

    .. attribute:: ave_locpot

        Mean of the minimum and maximmum (vacuum) locpot along c

    .. attribute:: sorted_sites

        List of sites from the slab sorted along the c direction

    .. attribute:: ave_bulk_p

        The average locpot of the slab region along the c direction
    """

    def __init__(self, poscar, locpot, outcar, shift=0):
        """
        Initializes the WorkFunctionAnalyzer class.

        Args:
            poscar (MSONable): Poscar vasp output object
            locpot (VolumetricData): Locpot vasp output object
            outcar (MSONable): Outcar vasp output object
            shift (float): Parameter to translate the slab (and
                therefore the vacuum) of the slab structure, thereby
                translating the plot along the x axis.
        """

        self.poscar = poscar
        self.locpot = locpot
        self.outcar = outcar

        # properties that can be shifted
        slab = poscar.structure.copy()
        slab.translate_sites([i for i, site in enumerate(slab)], [0, 0, shift])
        self.slab = slab
        self.sorted_sites = sorted(self.slab,
                                   key=lambda site: site.frac_coords[2])

        # Get the plot points between 0 and c
        # increments of the number of locpot points
        locpot_along_c = copy.copy(locpot.get_average_along_axis(2))
        self.along_c = np.linspace(0, 1, num=len(locpot_along_c))
        locpot_along_c_mid, locpot_end, locpot_start = [], [], []
        for i, s in enumerate(self.along_c):
            j = s + shift
            if j > 1:
                locpot_start.append(locpot_along_c[i])
            elif j < 0:
                locpot_end.append(locpot_along_c[i])
            else:
                locpot_along_c_mid.append(locpot_along_c[i])
        self.locpot_along_c = locpot_start + locpot_along_c_mid + locpot_end

        # get the average of the signal in the bulk-like region of the
        # slab, i.e. the average of the oscillating region. This gives
        # a rough appr. of the potential in the interior of the slab
        bulk_p = [p for i, p in enumerate(self.locpot_along_c) if \
                  self.sorted_sites[-1].frac_coords[2] > self.along_c[i] \
                  > self.sorted_sites[0].frac_coords[2]]
        self.ave_bulk_p = np.mean(bulk_p)

        # shift independent quantities
        self.efermi = outcar.efermi
        self.vacuum_locpot = max(self.locpot_along_c)
        # get the work function
        self.work_function = self.vacuum_locpot - self.efermi
        # for setting ylim and annotating
        self.ave_locpot = (self.vacuum_locpot - min(self.locpot_along_c)) / 2

    def get_locpot_along_slab_plot(self, label_energies=True, plt=None,
                                   label_fontsize=10):
        """
        Returns a plot of the local potential (eV) vs the
            position along the c axis of the slab model (Ang)

        Args:
            label_energies (bool): Whether to label relevant energy
                quantities such as the work function, Fermi energy,
                vacuum locpot, bulk-like locpot
            plt (plt): Matplotlib pylab object
            label_fontsize (float): Fontsize of labels

        Returns plt of the locpot vs c axis
        """

        plt = pretty_plot() if not plt else plt
        ax = list(plt.subplots())[1]
        ax.spines['top'].set_visible(False)

        # plot the raw locpot signal along c
        plt.plot(self.along_c, self.locpot_along_c, 'b--')

        # Get the local averaged signal of the locpot along c
        xg, yg = [], []
        for i, p in enumerate(self.locpot_along_c):
            # average signal is just the bulk-like potential when in the slab region
            if p < self.ave_bulk_p \
                    or self.sorted_sites[-1].frac_coords[2] >= self.along_c[i] \
                    >= self.sorted_sites[0].frac_coords[2]:
                yg.append(self.ave_bulk_p)
                xg.append(self.along_c[i])
            else:
                yg.append(p)
                xg.append(self.along_c[i])
        xg, yg = zip(*sorted(zip(xg, yg)))
        plt.plot(xg, yg, 'r', linewidth=2.5, zorder=-1)

        # make it look nice
        if label_energies:
            plt = self.get_labels(plt, label_fontsize=label_fontsize)
        plt.xlim([0, 1])
        plt.ylim([min(self.locpot_along_c),
                  self.vacuum_locpot + self.ave_locpot * 0.2])
        plt.xlabel(r"Fractional coordinates ($\hat{c}$)", fontsize=25)
        plt.xticks(fontsize=15, rotation=45)
        plt.ylabel(r"Energy (eV)", fontsize=25)
        plt.yticks(fontsize=15)

        return plt

    def get_labels(self, plt, label_fontsize=10):
        """
        Handles the optional labelling of the plot with relevant quantities
        Args:
            plt (plt): Plot of the locpot vs c axis
            label_fontsize (float): Fontsize of labels
        Returns Labelled plt
        """

        maxc = self.sorted_sites[-1].frac_coords[2]
        minc = self.sorted_sites[0].frac_coords[2]
        # determine whether most of the vacuum is to
        # the left or right for labelling purposes
        vleft = [i for i in self.along_c if i <= minc]
        vright = [i for i in self.along_c if i >= maxc]
        if max(vleft) - min(vleft) > max(vright) - min(vright):
            label_in_vac = (max(vleft) - min(vleft)) / 2
        else:
            label_in_vac = (max(vright) - min(vright)) / 2 + maxc

        # label the vacuum locpot
        label_in_bulk = maxc - (maxc - minc) / 2
        plt.plot([0, 1], [self.vacuum_locpot] * 2, 'b--', zorder=-5,
                 linewidth=1)
        xy = [label_in_bulk, self.vacuum_locpot + self.ave_locpot * 0.05]
        plt.annotate(r"$V_{vac}=%.2f$" % (self.vacuum_locpot), xy=xy,
                     xytext=xy, color='b', fontsize=label_fontsize)

        # label the fermi energy
        plt.plot([0, 1], [self.efermi] * 2, 'g--',
                 zorder=-5, linewidth=3)
        xy = [label_in_bulk, self.efermi + self.ave_locpot * 0.05]
        plt.annotate(r"$E_F=%.2f$" % (self.efermi), xytext=xy,
                     xy=xy, fontsize=label_fontsize, color='g')

        # label the bulk-like locpot
        plt.plot([0, 1], [self.ave_bulk_p] * 2, 'r--', linewidth=1., zorder=-1)
        print(label_in_vac)
        xy = [label_in_vac, self.ave_bulk_p + self.ave_locpot * 0.05]
        plt.annotate(r"$V^{interior}_{slab}=%.2f$" % (self.ave_bulk_p),
                     xy=xy, xytext=xy, color='r', fontsize=label_fontsize)

        # label the work function as a barrier
        plt.plot([label_in_vac] * 2, [self.efermi, self.vacuum_locpot],
                 'k--', zorder=-5, linewidth=2)
        xy = [label_in_vac, self.efermi + self.ave_locpot * 0.05]
        plt.annotate(r"$\Phi=%.2f$" % (self.work_function),
                     xy=xy, xytext=xy, fontsize=label_fontsize)

        return plt

    @staticmethod
    def from_files(poscar_filename, locpot_filename, outcar_filename, shift=0):
        poscar = Poscar.from_file(poscar_filename)
        locpot = Locpot.from_file(locpot_filename)
        outcar = Outcar(outcar_filename)
        return WorkFunctionAnalyzer(poscar, locpot, outcar, shift=shift)

    def as_dict(self):
        """

        Returns:

        """

        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "poscar": self.poscar.as_dict(),
             "locpot": self.locpot.as_dict(),
             "outcar": self.outcar.as_dict()}

        return d

    @classmethod
    def from_dict(cls, d):
        """

        Args:
            d:

        Returns:

        """

        return cls.__init__(d["poscar"], d["locpot"], d["outcar"])


class DielTensor(MSONable):
    """
    Class that represents the energy-dependent dielectric tensor of a solid
    state material.

    """

    def __init__(self, energies, dielectric_tensor):
        """
        Initializes a DielTensor instance from the dielectric data.

        Args:
            energies (numpy.array): (N,) array with the energy grid in eV.
            dielectric_tensor (numpy.array): (N, 3, 3) array with the dielectric
                tensor.

        """
        self._energies = energies
        self._dielectric_tensor = dielectric_tensor

    def check_dielectric_data(self):
        """
        Function that performs some tests on the dielectric data, to make sure
        input satisfies some constrains based on what we know about the dielectric
        tensor.

        Returns:
            None

        """
        pass  # TODO

    @property
    def energies(self):
        """
        Energy grid for which the dielectric tensor is defined in the original data.

        Returns:
            numpy.array: (N,) shaped array with the energies of the grid in eV.

        """
        return self._energies

    @property
    def dielectric_tensor(self):
        """
        Dielectric tensor of the material, calculated for each energy in the energy
        grid.

        Returns:
            numpy.array: (N, 3, 3) shaped array, where N corresponds to the number
                of energy grid points, and 3 to the different directions x,y,z.
        """
        return self._dielectric_tensor

    @property
    def dielectric_function(self):
        """
        The averaged dielectric function, derived from the tensor components by
        averaging the diagonal elements.

        Returns:
            np.array: (N,) shaped array with the dielectric function.
        """
        return np.array([np.mean(tensor.diagonal())
                         for tensor in self.dielectric_tensor])

    @property
    def absorption_coefficient(self):
        """
        Calculate the optical absorption coefficient from the dielectric data.
        For now the script only calculates the averaged absorption coefficient,
        i.e. by first averaging the diagonal elements and then using this
        dielectric function to calculate the absorption coefficient.

        Notes:
            The absorption coefficient is calculated as
            .. math:: \\alpha = \\frac{2 E}{ \hbar c} k(E)
            with $k(E)$ the imaginary part of the square root of the dielectric
            function

        Returns:
            np.array: (N,) shaped array with the energy (eV) dependent absorption
                coefficient in m^{-1}, where the energies correspond to self.energies.
        """

        energy = self.energies
        ext_coeff = np.array([cmath.sqrt(v).imag for v in self.dielectric_function])

        return 2.0 * energy * ext_coeff / (
                constants.hbar / constants.e * constants.c)

    def add_intraband_dieltensor(self, plasma_frequency, damping=0.1):
        """
        Add intraband component of the dielectric tensor based on the Drude model.

        Args:
            plasma_frequency:
            damping:

        Returns:

        """
        drude_diel = self.from_drude_model(plasma_frequency=plasma_frequency,
                                           energies=self.energies, damping=damping)

        self._dielectric_tensor += drude_diel.dielectric_tensor - 1

    def get_absorptivity(self, thickness, method="beer-lambert"):
        """
        Calculate the absorptivity for an absorber layer with a specified thickness
        and cell construction.

        Args:
            thickness (float): Thickness of the absorber layer, expressed in meters.
            method (str): Method for calculating the absorptivity.

        Returns:
            np.array: (N,) shaped array with the energy (eV) dependent absorptivity,
                where the energies correspond to self.energies.

        """
        if method == "beer-lambert":
            return 1.0 - np.exp(-2.0 * self.absorption_coefficient * thickness)
        else:
            raise NotImplementedError("Unrecognized method for calculating the "
                                      "absorptivity.")

    def get_loss_function(self, surface=False):
        """
        Calculate the loss function based on the averaged dielectric function of
        the dielectric tensor.

        Args:
            surface:

        Returns:

        """

        er = self.dielectric_function.real
        ei = self.dielectric_function.imag

        if surface:
            loss_function = ei / ((er + 1) ** 2 + ei ** 2)
        else:
            loss_function = ei / (er ** 2 + ei ** 2)

        loss_function[np.isnan(loss_function)] = 0

        return loss_function

    def plot(self, part="diel", variable_range=None, diel_range=None):
        """
        Plot the real and/or imaginary part of the dielectric function.

        Args:
            part (str): Which part of the dielectric function to plot, i.e. either
                "real", "imag" or "all".
            variable_range (tuple): Lower and upper limits of the variable which
                the requested function is plotted against.
            diel_range (tuple): Lower and upper limits of the range of the requested
                function in the plotted figure.

        Returns:
            None

        """
        if part == "diel":
            f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

            ax1.plot(self.energies, self.dielectric_function.real)
            ax2.plot(self.energies, self.dielectric_function.imag)
            if variable_range:
                ax1.set(xlim=variable_range)
                ax2.set(xlim=variable_range)
            if diel_range:
                ax1.set(ylim=diel_range)
                ax2.set(ylim=diel_range)
            ax1.set(ylabel=r"$\varepsilon_1$")
            ax2.set(xlabel="Energy (eV)", ylabel=r"$\varepsilon_2$")
            f.subplots_adjust(hspace=0.1)
            plt.show()

        elif part == "real":

            plt.plot(self.energies, self.dielectric_function.real)
            plt.xlabel("Energy (eV)")
            if variable_range:
                plt.xlim(variable_range)
            if diel_range:
                plt.ylim(diel_range)
            plt.ylabel(r"$\varepsilon_1$")
            plt.show()

        elif part == "imag":

            plt.plot(self.energies, self.dielectric_function.imag)
            plt.xlabel("Energy (eV)")
            if variable_range:
                plt.xlim(variable_range)
            if diel_range:
                plt.ylim(diel_range)
            plt.ylabel(r"$\varepsilon_2$")
            plt.show()

        elif part == "abs_coeff":

            plt.plot(self.energies, self.absorption_coefficient)
            plt.xlabel("Energy (eV)")
            if variable_range:
                plt.xlim(variable_range)
            if diel_range:
                plt.ylim(diel_range)
            plt.ylabel(r"$\alpha(E)$")
            plt.yscale("log")
            plt.show()

    def as_dict(self):
        """
        Note: stores the real and imaginary part of the dielectric tensor
        separately, due to issues with JSON serializing complex numbers.

        Returns:
            dict: Dictionary representation of the DielTensor instance.
        """
        d = dict()
        d["energies"] = MontyEncoder().default(self.energies)
        d["real_diel"] = MontyEncoder().default(self.dielectric_tensor.real)
        d["imag_diel"] = MontyEncoder().default(self.dielectric_tensor.imag)
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Initializes a DielTensor object from a dictionary.

        Args:
            d (dict): Dictionary from which the DielTensor should be initialized.

        Returns:
            DielTensor

        """
        energies = MontyDecoder().process_decoded(d["energies"])
        real_diel = MontyDecoder().process_decoded(d["real_diel"])
        imag_diel = MontyDecoder().process_decoded(d["imag_diel"])
        return cls(energies, real_diel + 1j * imag_diel)

    def to(self, filename):
        """
        Write the DielTensor to a JSON file.

        Args:
            filename (str): Path to the file in which the DielTensor should
                be written.

        Returns:
            None

        """
        with zopen(filename, "w") as f:
            f.write(self.to_json())

    @classmethod
    def from_file(cls, filename, fmt=None):
        """
        Initialize a DielTensor instance from a file.

        Args:
            filename (str): Path to file from which the dielectric data will be
                loaded. Can (so far) either be a vasprun.xml, OUTCAR or json file.
            fmt (str): Format of the file that contains the dielectric function
                data. Is optional, as the method can also figure out the format
                based on the filename.

        Returns:
            DielTensor: Dielectric tensor object from the dielectric data.

        """
        # Vasprun format: dielectric data is length 3 tuple
        if fmt == "vasprun" or filename.endswith(".xml"):
            dielectric_data = Vasprun(filename, parse_potcar_file=False).dielectric

            energies = np.array(dielectric_data[0])
            dielectric_tensor = np.array(
                [to_matrix(*real_data) + 1j * to_matrix(*imag_data)
                 for real_data, imag_data in zip(dielectric_data[1],
                                                 dielectric_data[2])]
            )
            return cls(energies, dielectric_tensor)

        # OUTCAR format: dielectric data is length 2 tuple
        elif fmt == "outcar" or fnmatch(filename, "*OUTCAR*"):
            outcar = Outcar(filename)
            outcar.read_freq_dielectric()
            return cls(outcar.frequencies, outcar.dielectric_tensor_function)

        # JSON format
        if fmt == "json" or filename.endswith(".json"):
            with zopen(filename, "r") as f:
                return cls.from_dict(json.loads(f.read()))

        else:
            raise IOError("Format of file not recognized. Note: Currently "
                          "only vasprun.xml and OUTCAR files are supported.")

    @classmethod
    def from_drude_model(cls, plasma_frequency, energies, damping=0.05):
        """
        Initialize a DielTensor object based on the Drude model for metals.

        Returns:

        """

        try:
            if not plasma_frequency.shape == (3, 3):
                raise ValueError("Plasma frequency array does not have right shape!")
        except AttributeError:
            if isinstance(plasma_frequency, float):
                plasma_frequency *= np.eye(3)
            else:
                raise TypeError("The plasma frequency must be expressed either as "
                                "a float or a numpy array of shape (3, 3).")

        dieltensor = np.array(
            [1 - omega ** 2 / (energies ** 2 + 1j * energies * damping)
             for omega in plasma_frequency.reshape(9)]
        ).reshape((3, 3, len(energies)))
        dieltensor = dieltensor.swapaxes(0, 2)

        return cls(energies, dieltensor)


# Utility method

def to_matrix(xx, yy, zz, xy, yz, xz):
    """
    Convert a list of matrix components to a symmetric 3x3 matrix.
    Inputs should be in the order xx, yy, zz, xy, yz, xz.

    Args:
        xx (float): xx component of the matrix.
        yy (float): yy component of the matrix.
        zz (float): zz component of the matrix.
        xy (float): xy component of the matrix.
        yz (float): yz component of the matrix.
        xz (float): xz component of the matrix.

    Returns:
        (np.array): The matrix, as a 3x3 numpy array.

    """
    matrix = np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])
    return matrix


def fix_slab_bulk(poscar, thickness, method="layers", part="center"):
    """
    Fix atoms of a slab to represent the bulk of the material. Which atoms are
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

    Returns:
        selective dynamics (Nx3 array): bool values for selective dynamics,
            where N is number of sites. Defaults to None.

    """
    # TODO : This class has been moved elsewhere. Maybe I should remove it here,
    #  or it could be useful for others...

    if method == "layers":

        atomic_layers = find_atomic_layers(poscar.structure)

        if part == "center":

            # Even number of layers
            if len(atomic_layers) % 2 == 0:

                # Check if the user requested an odd number of layers for the
                # fixed part of the slab
                if thickness % 2 == 1:
                    print("Found an even number of layers, but the user " +
                          "requested an odd number of fixed layers. Adding "
                          "one layer to the fixed part of the slab.")
                    thickness += 1

            # Odd number of layers
            if len(atomic_layers) % 2 == 1:

                # Check if the user requested an even number of layers for the
                # fixed part of the slab
                if thickness % 2 == 0:
                    print("Found an odd number of layers, but the user " +
                          "requested an even number of fixed layers. Adding "
                          "one layer to the fixed part of the slab.")
                    thickness += 1

            # Calculate the number of layers to optimize on each site
            n_optimize_layers = int((len(atomic_layers) - thickness) / 2)

            if n_optimize_layers < 5:
                print("WARNING: Less than 5 layers are optimized on each "
                      "side of the slab.")

            # Take the fixed layers from the atomic layers of the slab
            fixed_layers = atomic_layers[n_optimize_layers: n_optimize_layers +
                                                            thickness]

        else:
            raise NotImplementedError("Requested part is not implemented " +
                                      "(yet).")
            # TODO Implement oneside

        fixed_sites = [site for layer in fixed_layers for site in layer]

        # Set up the selective dynamics property

        selective_dynamics = []

        for site in poscar.structure.sites:
            if site in fixed_sites:
                selective_dynamics.append([False, False, False])
            else:
                selective_dynamics.append([True, True, True])

        return selective_dynamics

    else:
        raise NotImplementedError("Requested method is not implemented (yet).")
        # TODO Implement angstrom


def write_all_slab_terminations(structure, miller_indices, min_slab_size,
                                min_vacuum_size):
    """
    Writes the POSCAR files for all slab terminations.

    Args:
        structure:
        miller_indices:
        min_slab_size:
        min_vacuum_size:

    Returns:
        None
    """

    slab_gen = SlabGenerator(structure, miller_indices, min_slab_size,
                             min_vacuum_size)
    slabs = slab_gen.get_slabs()

    letter_counter = 0

    for slab in slabs:
        slab.sort(reverse=True)

        slab_letter = string.ascii_lowercase[letter_counter]
        letter_counter += 1

        filename = "".join([str(number) for number in miller_indices]) + "_" \
                   + slab_letter + "_" + str(min_slab_size) + "l_POSCAR.vasp"

        slab.to(fmt='vasp', filename=filename)

        print("Written slab structure to " + filename)
        print("Slab is symmetric = " + str(slab.is_symmetric()))
        print("Slab is polar = " + str(slab.is_polar()))


def find_suitable_kpar(structure, kpoints, max_kpar=30):
    """

    :param structure:
    :param kpoints:
    :return:
    """

    spg = SpacegroupAnalyzer(structure)

    kpar = len(spg.get_ir_reciprocal_mesh(kpoints.kpts))
    divisors = generate_divisors(kpar)

    while kpar > max_kpar:
        kpar = next(divisors)

    return kpar


def find_irr_kpoints(structure, kpoints):
    """

    :param structure:
    :param kpoints:
    :return:
    """

    spg = SpacegroupAnalyzer(structure, symprec=1e-5)

    return len(spg.get_ir_reciprocal_mesh(kpoints.kpts))


def generate_divisors(number):
    """

    Args:
        number:

    Returns:

    """
    divisor = int(number / 2)

    while divisor != 0:

        while number % divisor != 0:
            divisor -= 1

        yield divisor
        divisor -= 1


# Stolen from https://goshippo.com/blog/measure-real-size-any-python-object/
# Used to calculate the actual total size of a python object


def get_size(obj, seen=None):
    """Recursively finds size of objects"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj,
                                                     (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size
