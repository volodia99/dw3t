from dataclasses import dataclass
import os
from typing import Any
from pathlib import Path
import urllib.request

import numpy as np
import astropy.units as u
import astropy.constants as uc

import inifix
import dsharp_opac as do

from nonos._geometry import axes_from_geometry, Geometry
from dw3t._typing import FArray1D, FArrayND
from dw3t._parsing import is_set

@dataclass(kw_only=True, slots=True, frozen=True)
class Grid:
    x1: FArray1D[u.Quantity]
    x2: FArray1D[u.Quantity]
    x3: FArray1D[u.Quantity]

    @property
    def x1c(self) -> FArray1D[u.Quantity]:
        return 0.5*(self.x1[1:]+self.x1[:-1])

    @property
    def x2c(self) -> FArray1D[u.Quantity]:
        return 0.5*(self.x2[1:]+self.x2[:-1])

    @property
    def x3c(self) -> FArray1D[u.Quantity]:
        return 0.5*(self.x3[1:]+self.x3[:-1])

    @property
    def shape(self):
        return (self.x1c.shape[0], self.x2c.shape[0], self.x3c.shape[0])

@dataclass(kw_only=True, slots=True, frozen=True)
class Dust:
    rho: FArrayND[u.Quantity]
    size: u.Quantity

@dataclass(kw_only=True, slots=True, frozen=True)
class Gas:
    rho: FArrayND[u.Quantity]
    v1: FArrayND[u.Quantity]
    v2: FArrayND[u.Quantity]
    v3: FArrayND[u.Quantity]

    @property
    def velocity(self):
        return np.stack((self.v1, self.v2, self.v3), axis=0)

#TODO: Opacity to be improved
@dataclass(kw_only=True, slots=True)
class Opacity:
    mix:str
    rho:float|None=None

    def __post_init__(self):
        if self.mix.endswith(".lnk"):
            self.mix = do.diel_from_lnk_file(self.mix)
            if not is_set(self.rho):
                raise ValueError(
                    f"Internal density of the mix has to be defined. Please provide 'rho' in dust.opacity."
                )
            self.mix.rho = self.rho
        else:
            if is_set(self.rho):
                print("WARNING: unused 'rho' when using dsharp_opac mix.")

@dataclass(kw_only=True, slots=True)
class Model:
    grid:Grid|None=None
    dust:Dust|None=None
    gas:Gas|None=None
    unit_length_au:float
    unit_mass_msun:float
    component:str
    geometry:Geometry

    @property
    def dimension(self):
        return (3-self.grid.shape.count(1))

    @property
    def axes_from_geometry(self):
        return axes_from_geometry(self.geometry)

    @property
    def _indices_of_reduced_axes(self):
        if self.axes_from_geometry.count(1) == 0:
            indices = None
        elif self.axes_from_geometry.count(1) == 1:
            indices = self.axes_from_geometry.index(1)
        elif self.axes_from_geometry.count(1) > 1:
            indices = [i for i, x in enumerate(self.axes_from_geometry) if x == 1]
        return(np.atleast_1d(indices))

    @property
    def reduced_axes(self):
        if None not in self._indices_of_reduced_axes:
            return tuple(np.array(self.axes_from_geometry[self._indices_of_reduced_axes]))
        else:
            return (None,)

    def write_files(
        self, 
        *, 
        directory:str,
        write_opacities:bool=False,
        opacity=None,
        smoothing:bool=False,
        simulation_files_only:bool=False,
        binary:bool=False,
        config:dict,
        ):
        """
        Adapted from dustpylib (https://dustpylib.readthedocs.io/en/latest/radmc3d.html) method 
        Function writes all required ``RADMC-3D`` input files.

        Parameters
        ----------
        directory : str
            Data directory in which the files are written.
        """
        Path(directory).mkdir(parents=True, exist_ok=False)
        if simulation_files_only:
            print(f"WARNING: for now writing only the input files related to simulated outputs, "\
            "refer to the RADMC3D documentation to write the mandatory files related to RADMC3D.")
        elif not simulation_files_only:
            self._write_radmc3d_inp(directory=directory, config=config["radmc3d"])
            self._write_stars_inp(directory=directory, config=config)
            self._write_wavelength_micron_inp(directory=directory, config=config["wavelength_micron"])
            if "gas" in self.component:
                self._write_lines_inp(directory=directory, config=config["gas"])
                self._write_molecule_inp(directory=directory, config=config["gas"])
        self._write_amr_grid_inp(directory=directory)
        if "dust" in self.component:
            self._write_dust_density_inp(directory=directory)
            if write_opacities:
                self.write_opacity_files(
                    directory=directory,
                    opacity=opacity.mix,
                    smoothing=smoothing,
                    config=config,
                )
            self._write_metadata(directory=directory)
        if "gas" in self.component:
            self._write_numberdens_inp(directory=directory, config=config)
            self._write_gas_velocity_inp(directory=directory)

    def _write_radmc3d_inp(self, *, directory:str, config:dict):
        """
        Adapted from dustpylib (https://dustpylib.readthedocs.io/en/latest/radmc3d.html) method 
        Function writes the 'radmc3d.inp' input file.

        Parameters
        ----------
        directory : str
            Data directory in which the files are written.
        """

        filename = "radmc3d.inp"
        path = os.path.join(directory, filename)

        print(f"INFO: Writing {path}.....", end="")
        with open(path, "w") as f:
            for key in config:
                f.write(f"{key} = {config[key]}\n")
        print("done.")

    def _write_stars_inp(self, *, directory:str, config:dict):
        """
        Adapted from dustpylib (https://dustpylib.readthedocs.io/en/latest/radmc3d.html) method 
        Function writes the 'stars.inp' input file.

        Parameters
        ----------
        directory : str
            Data directory in which the files are written.
        """

        filename = "stars.inp"
        path = os.path.join(directory, filename)

        # Set up a wavelength grid (in micron) upon which we want to compute the opacities
        config_wavelength_micron = config["wavelength_micron"]
        lam_grid = np.linspace(config_wavelength_micron["min"],config_wavelength_micron["max"],config_wavelength_micron["N"])*u.micron
        lam_grid = lam_grid.value
        R_star = (config["stars"]["R_star"]*u.R_sun).to(u.cm).value
        M_star = (config["stars"]["M_star"]*u.M_sun).to(u.g).value
        T_star = (config["stars"]["T_star"]*u.K).value

        print(f"INFO: Writing {path}.....", end="")
        with open(path, "w") as f:
            f.write("2\n")
            f.write(f"1 {lam_grid.shape[0]:d}\n")
            f.write(
                f"{R_star:.6e} "\
                f"{M_star:.6e} "\
                f"{0:.6e} "\
                f"{0:.6e} "\
                f"{0:.6e}\n"
                )
            for lam in lam_grid:
                f.write(f"{lam:.6e}\n")
            f.write(f"-{T_star:.6e}\n")
        print("done.")

    def _write_wavelength_micron_inp(self, *, directory:str, config:dict):
        """
        Adapted from dustpylib (https://dustpylib.readthedocs.io/en/latest/radmc3d.html) method 
        Function writes the 'wavelength_micron.inp' input file.

        Parameters
        ----------
        directory : str
            Data directory in which the files are written.
        """

        filename = "wavelength_micron.inp"
        path = os.path.join(directory, filename)

        # Set up a wavelength grid (in micron) upon which we want to compute the opacities
        lam_grid = np.linspace(config["min"],config["max"],config["N"])*u.micron
        lam_grid = lam_grid.value

        print(f"INFO: Writing {path}.....", end="")
        with open(path, "w") as f:
            f.write(f"{lam_grid.shape[0]:d}\n")
            for lam in lam_grid:
                f.write(f"{lam:.6e}\n")
        print("done.")

    def _write_lines_inp(self, *, directory:str, config:dict):
        """
        Function writes the 'lines.inp' input file.

        Parameters
        ----------
        directory : str
            Data directory in which the files are written.
        """
        filename = "lines.inp"
        path = os.path.join(directory, filename)

        species = config["species"]
        #TODO: be more flexible for NON-LTE calculations?
        non_lte = 0

        print(f"INFO: Writing {path}.....", end="")
        with open(path, "w") as f:
            f.write("2\n")
            f.write("1\n")
            f.write(
                f"{species} "\
                "leiden "\
                f"{0:d} "\
                f"{0:d} "\
                f"{non_lte:d}\n"
                )
        print("done.")

    def _write_molecule_inp(self, *, directory:str, config:dict):
        """
        Function writes the 'molecule_*.inp' input file.

        Parameters
        ----------
        directory : str
            Data directory in which the files are written.
        """
        species = config["species"]
        filename = f"molecule_{species}.inp"
        path = os.path.join(directory, filename)

        # print(f"INFO: Writing {path}.....", end="")
        urllib.request.urlretrieve(
            f"https://home.strw.leidenuniv.nl/~moldata/datafiles/{species}.dat", 
            path,
        )

    def _write_metadata(self, *, directory:str|None=None):
        """
        Adapted from dustpylib (https://dustpylib.readthedocs.io/en/latest/radmc3d.html) method 
        Function writes the 'metadata.npz' file.
        It is not required for ``RADMC-3D`` model, but it contains
        additional data such as particle size bins used in the setup.

        Parameters
        ----------
        directory : str
            Data directory in which the files are written.
        """

        filename = "metadata.dw3t.npz"
        if directory is None:
            raise ValueError("directory for RT calculation should be specified.")
        path = os.path.join(directory, filename)

        print(f"INFO: Writing {path}.....", end="")

        np.savez(
            path,
            dust_size_array=self.dust.size.value,
            dust_size_unit=str(self.dust.size.unit),
        )
        print("done.")

    def write_opacity_files(
            self, *, directory:str, opacity, smoothing:bool=False, config:dict):
        """
        Adapted from dustpylib (https://dustpylib.readthedocs.io/en/latest/radmc3d.html) method 
        Function writes the required opacity files.

        Parameters
        ----------
        directory : str
            Data directory in which the files are written.
        opacity : str or instance of `dharp_opac.diel_const`
            Opacity model to be used. Either 'birnstiel2018', 'ricci2010', or instance of do.diel_const.
            Could be read from a .lnk file using do.diel_from_lnk_file class, inherited from do.diel_const.
            See https://github.com/birnstiel/dsharp_opac for more details 
        smoothing : bool, optional, default: False
            Smooth the opacities by averaging over multiple particle sizes.
            This slows down the computation.
        """
        self._write_dustopac_inp(directory=directory)
        self._write_dustkapscatmat_inp(
            directory=directory,
            opacity=opacity,
            smoothing=smoothing,
            config=config["wavelength_micron"],
        )

    def _write_dustopac_inp(self, *, directory:str):
        """
        Adapted from dustpylib (https://dustpylib.readthedocs.io/en/latest/radmc3d.html) method 
        Function writes the 'dustopac.inp' input file.

        Parameters
        ----------
        directory : str
            Data directory in which the files are written.
        """

        filename = "dustopac.inp"
        path = os.path.join(directory, filename)

        print(f"INFO: Writing {path}.....", end="")
        Nspec = self.dust.size.shape[0]
        # computes the 'magnitude' of Nspec
        mag = int(np.ceil(np.log10(Nspec)))
        with open(path, "w") as f:
            f.write("2\n")
            f.write(f"{Nspec}\n")
            for idust in range(Nspec):
                f.write("--------------------\n")
                f.write("10\n")
                f.write("0\n")
                f.write(f"{str(idust).zfill(mag)}\n")
        print("done.")

    def _write_dustkapscatmat_inp(
        self, 
        *, 
        directory:str,
        opacity,
        smoothing:bool=False,
        config:dict,
        ):
        """
        Adapted from dustpylib (https://dustpylib.readthedocs.io/en/latest/radmc3d.html) method 
        Function writes the 'dustkapscatmat_*.inp' input files.

        Parameters
        ----------
        directory : str
            Data directory in which the files are written.
        opacity : str or instance of `dharp_opac.diel_const`
            Opacity model to be used. Either 'birnstiel2018' or 'ricci2010' or instance of do.diel_const.
            Could be read from a .lnk file using do.diel_from_lnk_file class, inherited from do.diel_const.
        smoothing : bool, optional, default: False
            Smooth the opacities by averaging over multiple particle sizes.
            This slows down the computation.
        """
        #TODO: be more flexible?
        Nangle = 181
        # Set up a wavelength grid (in micron) upon which we want to compute the opacities
        lam_grid = np.linspace(config["min"],config["max"],config["N"])*u.micron
        Nlam = lam_grid.shape[0]

        dustsize = self.dust.size.to(u.cm).value
        Nspec = dustsize.shape[0]
        mag = int(np.ceil(np.log10(Nspec)))

        print()
        print("INFO: Computing opacities...")
        print("INFO: Using dsharp_opac. Please cite Birnstiel et al. (2018).")

        # Selecting the opacity model
        if opacity == "birnstiel2018":
            print("INFO: Using DSHARP mix. Please cite Birnstiel et al. (2018).")
            mix, rho_s = do.get_dsharp_mix()
        elif opacity == "ricci2010":
            print("INFO: Using Ricci mix. Please cite Ricci et al. (2010).")
            mix, rho_s = do.get_ricci_mix(lmax=lam_grid[-1].to(u.cm).value,
                                          extrapol=True)
        elif isinstance(opacity, do.diel_const):
            mix = opacity
            if hasattr(opacity, 'rho') and opacity.rho is not None:
                rho_s = opacity.rho
            else:
                raise ValueError(
                    'opacity needs to have the attribute rho (material density) set')
        else:
            raise RuntimeError(f"Unknown opacity '{opacity}'")

        # When smoothing is True the opacities are computed on a finer grid (x4)
        if smoothing:
            amin = dustsize.min()
            amax = dustsize.max()
            #TODO: be more flexible?
            Na = 4*dustsize.shape[0]
            a_opac = np.geomspace(amin, amax, Na)
        else:
            a_opac = dustsize

        # Computing the opacities
        opac_dict = do.get_opacities(
            a_opac, lam_grid.to(u.cm).value,
            rho_s, mix,
            extrapolate_large_grains=True,
            n_angle=int((Nangle-1)/2+1)
        )

        # When smoothing is True several size bins are averaged into
        # the actual RADMC-3D model bins
        if smoothing:
            k_abs = np.empty((dustsize.shape[0], lam_grid.shape[0]))
            k_sca = np.empty((dustsize.shape[0], lam_grid.shape[0]))
            g_sca = np.empty((dustsize.shape[0], lam_grid.shape[0]))
            S1 = np.empty(
                (dustsize.shape[0], lam_grid.shape[0], Nangle),
                dtype=complex)
            S2 = np.empty(
                (dustsize.shape[0], lam_grid.shape[0], Nangle),
                dtype=complex)
            for i in range(dustsize.shape[0]):
                # This is equivalent to weighting the opacities with an MRN
                # size distribution of n(a) \protp a^{-3.5}
                sqrta = np.sqrt(a_opac)
                sqrta_sum = sqrta.sum()
                k_abs[i, :] = (sqrta[:, None] * opac_dict["k_abs"]
                               ).mean(0) / sqrta_sum
                k_sca[i, :] = (sqrta[:, None] * opac_dict["k_sca"]
                               ).mean(0) / sqrta_sum
                g_sca[i, :] = (sqrta[:, None] * opac_dict["g"]
                               ).mean(0) / sqrta_sum
                S1[i, ...] = (sqrta[:, None, None] * opac_dict["S1"]
                              ).mean(0) / sqrta_sum
                S2[i, ...] = (sqrta[:, None, None] * opac_dict["S2"]
                              ).mean(0) / sqrta_sum
            opac_dict["k_abs"] = k_abs
            opac_dict["k_sca"] = k_sca
            opac_dict["g"] = g_sca
            opac_dict["S1"] = S1
            opac_dict["S2"] = S2
            opac_dict["a"] = dustsize

        # Making sure that the scattering phase functions are
        # properly normalized.
        zscat, _, k_sca, g = do.chop_forward_scattering(opac_dict)
        opac_dict["k_sca"] = k_sca
        opac_dict["g"] = g
        opac_dict["zscat"] = zscat
        print()

        # Writing the files
        for ia in range(Nspec):
            filename = f"dustkapscatmat_{str(ia).zfill(mag)}.inp"
            path = os.path.join(directory, filename)
    
            print(f"INFO: Writing {path}.....", end="")
            with open(path, "w") as f:
                f.write("1\n")
                f.write(f"{Nlam:d}\n")
                f.write(f"{Nangle:d}\n")
                f.write("\n")
                for ilam in range(Nlam):
                    f.write(
                        f"{lam_grid[ilam].value:.6e} "\
                        f"{opac_dict["k_abs"][ia, ilam]:.6e} "\
                        f"{opac_dict["k_sca"][ia, ilam]:.6e} "\
                        f"{opac_dict["g"][ia, ilam]:.6e}\n"
                    )
                f.write("\n")
                for theta in opac_dict["theta"]:
                    f.write(f"{theta:.2f}\n")
                f.write("\n")
                for ilam in range(Nlam):
                    for iang in range(Nangle):
                        f.write(
                            f"{zscat[ia, ilam, iang, 0]:.6e} "\
                            f"{zscat[ia, ilam, iang, 1]:.6e} "\
                            f"{zscat[ia, ilam, iang, 2]:.6e} "\
                            f"{zscat[ia, ilam, iang, 3]:.6e} "\
                            f"{zscat[ia, ilam, iang, 4]:.6e} "\
                            f"{zscat[ia, ilam, iang, 5]:.6e}\n"
                        )
            print("done.")
        print()

    def _write_numberdens_inp(self, *, directory:str, config:dict):
        """
        Function writes the 'numberdens_*.inp' input file.

        Parameters
        ----------
        directory : str
            Data directory in which the files are written.
        """
        config_gas = config["gas"]
        species = config_gas["species"]
        abundance = config_gas["abundance"]
        if abundance["mode"] not in ("constant", "array"):
            raise ValueError(f"abundance.mode = {abundance["mode"]}. Should be 'constant' or 'array' from npz file.")
        filename = f"numberdens_{species}.binp"
        path = os.path.join(directory, filename)

        MUSTAR = config["stars"]["mu_star"]
        numberrho_H2 = self.gas.rho/(MUSTAR*uc.m_p.to(u.g))
        if numberrho_H2.unit!=1/u.cm**3:
            raise ValueError(f"gas rho field unit should be 1/cm^3, not {numberrho_H2.unit}.")
        if abundance["mode"]=="constant":
            numberrho = numberrho_H2*abundance["value"]
        elif abundance["mode"]=="array":
            raise NotImplementedError("mode='array' not yet implemented.")
        print(f"INFO: Writing {path}.....", end="")
        with open(path, "wb") as f:
            header = np.array([1, 8, numberrho.flatten().shape[0]], dtype=int)
            header.tofile(f)
            numberrho.ravel(order="F").value.tofile(f)
            # for ngrho in numberrho.ravel(order="F"):
            #     f.write(f"{ngrho.value:.6e}\n".encode("utf-8"))
        print("done.")

    def _write_gas_velocity_inp(self, *, directory:str):
        """
        Function writes the 'gas_velocity.inp' input file.

        Parameters
        ----------
        directory : str
            Data directory in which the files are written.
        """
        filename = "gas_velocity.binp"
        if len(set({self.gas.v1.unit, self.gas.v2.unit, self.gas.v3.unit}))!=1 and self.gas.v1.unit!=u.cm/u.s:
            raise ValueError(f"gas velocity fields unit should be cm/s, not {self.gas.v1.unit}.")
        path = os.path.join(directory, filename)

        print(f"INFO: Writing {path}.....", end="")
        with open(path, "wb") as f:
            header = np.array([1, 8, self.gas.v1.flatten().shape[0]], dtype=int)
            header.tofile(f)
            self.gas.velocity.ravel(order="F").value.tofile(f)
            # for gasv1, gasv2, gasv3 in zip(self.gas.v1.ravel(order="F"), self.gas.v2.ravel(order="F"), self.gas.v3.ravel(order="F")):
            #     f.write((
            #         f"{gasv1.value:.6e} "\
            #         f"{gasv2.value:.6e} "\
            #         f"{gasv3.value:.6e}\n").encode("utf-8")
            #     )
        print("done.")

    def _write_amr_grid_inp(self, *, directory:str):
        """
        Adapted from dustpylib (https://dustpylib.readthedocs.io/en/latest/radmc3d.html) method 
        Function writes the 'amr_grid.inp' input file.

        Parameters
        ----------
        directory : str
            Data directory in which the files are written.
        """

        filename = "amr_grid.inp"
        if self.geometry!="spherical":
            raise NotImplementedError(f"Model geometry={self.geometry} should be 'spherical'.")
        if self.grid.x1.unit!=u.cm:
            raise ValueError(f"radial coordinate unit should be cm, not {self.grid.x1.unit}.")
        path = os.path.join(directory, filename)

        print(f"INFO: Writing {path}.....", end="")
        Nr = self.grid.x1c.shape[0]
        Ntheta = self.grid.x2c.shape[0]
        Nphi = self.grid.x3c.shape[0]
        with open(path, "w") as f:
            f.write("1\n")
            f.write("0\n")
            f.write("100\n")
            f.write("0\n")
            f.write(
                f"{(1 if Nr > 1 else 0):d} "\
                f"{1 if Ntheta > 1 else 0:d} "\
                f"{1 if Nphi > 1 else 0:d}\n"
            )
            f.write(f"{Nr:d} {Ntheta:d} {Nphi:d}\n")
            for x1 in self.grid.x1:
                f.write(f"{x1.value:.12e}\n")
            for x2 in self.grid.x2:
                f.write(f"{x2.value:.12e}\n")
            for x3 in self.grid.x3:
                f.write(f"{x3.value:.12e}\n")
        print("done.")

    def _write_dust_density_inp(self, *, directory:str):
        """
        Adapted from dustpylib (https://dustpylib.readthedocs.io/en/latest/radmc3d.html) method 
        Function writes the 'dust_density.inp' input file.

        Parameters
        ----------
        directory : str
            Data directory in which the files are written. 
        """
        filename = "dust_density.binp"
        if self.dust.rho.unit!=u.g/u.cm**3:
            raise ValueError(f"dust rho field unit should be g/cm^3, not {self.dust.rho.unit}.")
        path = os.path.join(directory, filename)

        print(f"INFO: Writing {path}.....", end="")
        with open(path, "wb") as f:
            header = np.array([1, 8, self.dust.rho[..., 0].flatten().shape[0], self.dust.size.shape[0]], dtype=int)
            header.tofile(f)
            self.dust.rho.ravel(order="F").value.tofile(f)
            # for dustrho in self.dust.rho.ravel(order="F"):
            #     f.write(f"{dustrho.value:.6e}\n".encode("utf-8"))
        print("done.")

def load_model(
    unit_length_au:float,
    unit_mass_msun:float,
    component:list,
    config:dict,
) -> "Model":
    if not (set(component) & set(["gas","dust"])):
        raise ValueError(
            f"{component=} should be 'gas', 'dust' or ['dust', 'gas']."
        )
    unit_length_au = unit_length_au * u.au
    unit_mass_msun = unit_mass_msun * u.M_sun

    model = Model(
        unit_length_au=unit_length_au,
        unit_mass_msun=unit_mass_msun,
        component=component,
        geometry=Geometry("spherical")
    )
    return model
