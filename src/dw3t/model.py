from dataclasses import dataclass, field
import os
from typing import Any
from pathlib import Path
import urllib.request
import sys


import numpy as np
import astropy.units as u
import astropy.constants as uc
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

import inifix
import dsharp_opac as do
import prodimopy.read as pread

from nonos._geometry import axes_from_geometry, Geometry
from dw3t._typing import FArray1D, FArrayND
from dw3t._parsing import is_set

if sys.version_info >= (3, 11):
    from typing import assert_never
else:
    from typing_extensions import assert_never

if sys.version_info >= (3, 13):
    from copy import replace
else:
    from dataclasses import replace


@dataclass(kw_only=True, slots=True, frozen=True)
class CellCenters3D:
    x1c: FArray1D[u.Quantity]
    x2c: FArray1D[u.Quantity]
    x3c: FArray1D[u.Quantity]
    geometry: Geometry

    def __post_init__(self):
        if self.x1c.shape!=self.x2c.shape!=self.x3c.shape:
            raise ValueError(
                f"{self.x1c.shape=}, {self.x2c.shape=} and {self.x3c.shape=} should be the same."
            )

    @property
    def shape(self):
        return self.x1c.shape

@dataclass(slots=True, frozen=True, kw_only=True)
class Array:
    cells:CellCenters3D
    data:FArrayND
    config:dict

    def __post_init__(self):
        if self.cells.shape!=self.data.shape:
            raise ValueError(
                f"{self.cells.shape=} should be equal to {self.data.shape=}"
            )
        try:
            unit = self.data.unit
        except AttributeError:
            raise AttributeError(
                f"data should have a unit."
            )

    def _cart2pol(self, *, x, z):
        """
        Conversion of the coordinate system to cartesian.
        """
        r = np.sqrt(x**2+z**2)
        theta = np.arctan2(x, z)
        return (r, theta)

    def interpolate_on_half_spherical_disk(self, output_cells:"CellCenters3D") -> "Array":
        if output_cells.geometry!="spherical":
            raise NotImplementedError(f"Grid geometry={output_cells.geometry} should be 'spherical'.")
        #TODO: add raise condition if interpolate_on_half_spherical_disk not applicable, ie not polar coordinate and 2D ?
        unit_length_au = self.config["simulation"]["unit_length_au"]

        x_3d = np.flip(self.cells.x1c, axis=1)
        z_3d = np.flip(self.cells.x2c, axis=1)
        phi_3d = np.flip(self.cells.x3c, axis=1)
        data_species = np.flip(self.data, axis=1)
        r_3d, theta_3d = self._cart2pol(
            x=x_3d,
            z=z_3d,
        )
    
        ntheta = output_cells.shape[1]
        output_r2d_half = output_cells.x1c[:,0:ntheta//2,0]
        output_theta2d_half = output_cells.x2c[:,0:ntheta//2,0]

        points = np.array([(r_3d[:,:,0].to(u.au)/(unit_length_au*u.au)).value.flatten(), theta_3d[:,:,0].value.flatten()]).T
        interpolated_data_species = griddata(points, data_species[:,:,0].flatten(), ((output_r2d_half.to(u.au)/(unit_length_au*u.au)).value, output_theta2d_half.value), method="linear")
        return replace(
            self,
            cells=CellCenters3D(
                x1c=np.atleast_3d(output_r2d_half),
                x2c=np.atleast_3d(output_theta2d_half),
                x3c=np.zeros_like(np.atleast_3d(output_r2d_half).value)*u.radian,
                geometry=output_cells.geometry,
            ),
            data=np.atleast_3d(interpolated_data_species)*data_species.unit,
        )

    def phi_expansion(self, output_cells:"CellCenters3D") -> "Array":
        if output_cells.geometry!="spherical":
            raise NotImplementedError(f"Grid geometry={output_cells.geometry} should be 'spherical'.")
        if (phi_shape:=self.data.shape[2])!=1:
            raise ValueError(
                f"{phi_shape=}. It should be 1 to perform a phi_expansion operation."
            )
        if self.data.shape!=self.cells.x3c.shape:
            raise ValueError(
                f"There is a discrepancy between the data shape: {self.data.shape}, and cells.x3c: {self.cells.x3c.shape}."
            )

        processed_data = np.broadcast_to(
            self.data, 
            output_cells.shape,
        )*self.data.unit
        processed_x1c = np.broadcast_to(
            self.cells.x1c, 
            output_cells.shape,
        )*self.cells.x1c.unit
        processed_x2c = np.broadcast_to(
            self.cells.x2c, 
            output_cells.shape,
        )*self.cells.x2c.unit
        processed_cells = CellCenters3D(
            x1c=processed_x1c,
            x2c=processed_x2c,
            x3c=output_cells.x3c,
            geometry=output_cells.geometry,
        )
        return replace(
            self,
            cells=processed_cells,
            data=processed_data,
        )

    def mirror_symmetry(self, output_cells:"CellCenters3D") -> "Array":
        if output_cells.geometry!="spherical":
            raise NotImplementedError(f"Grid geometry={output_cells.geometry} should be 'spherical'.")
        #TODO: add raise condition if mirror_symmetry not applicable
        processed_data = np.concatenate(
            [
                self.data,
                np.flip(self.data, axis=1),
            ],
            axis=1,
        )
        ntheta = output_cells.shape[1]
        processed_x1c = np.concatenate(
            [
                self.cells.x1c,
                np.atleast_3d(output_cells.x1c[:,ntheta//2:ntheta,0])
            ],
            axis=1,
        )
        processed_x3c = np.concatenate(
            [
                self.cells.x3c,
                np.atleast_3d(output_cells.x3c[:,ntheta//2:ntheta,0])
            ],
            axis=1,
        )

        processed_cells = CellCenters3D(
            x1c=processed_x1c,
            x2c=np.atleast_3d(output_cells.x2c[:,:,0]),
            x3c=processed_x3c,
            geometry=output_cells.geometry,
        )
        if processed_data.shape!=processed_cells.shape:
            raise ValueError(
                f"There is a discrepancy between the data shape: {processed_data.shape}, and the grid shape: {processed_cells.shape}."
            )
        return replace(
            self,
            cells=processed_cells,
            data=processed_data,
        )
    
    def remove_nan(self) -> "Array":
        processed_data = np.nan_to_num(self.data, nan=np.nanmin(self.data))
        return replace(
            self,
            data=processed_data,
        )

    def smooth_gaussian_filter(self) -> "Array":
        processed_data = gaussian_filter(self.data, sigma=1)*self.data.unit
        return replace(
            self,
            data=processed_data,
        )

    def from_prodimo_to(self, output_cells:"CellCenters3D") -> "Array":
        #TODO: all operations in one --> lot of work that can't be used anyway
        interpolated_array = self.interpolate_on_half_spherical_disk(output_cells=output_cells)
        mirrored_array = interpolated_array.mirror_symmetry(output_cells=output_cells)
        phi_expanded_array = mirrored_array.phi_expansion(output_cells=output_cells)
        no_nan_array = phi_expanded_array.remove_nan()
        smoothed_array = no_nan_array.smooth_gaussian_filter()
        array = smoothed_array
        return replace(
            self,
            cells=output_cells,
            data=array.data,
        )

@dataclass(kw_only=True, slots=True, frozen=True)
class Grid:
    x1: FArray1D[u.Quantity]
    x2: FArray1D[u.Quantity]
    x3: FArray1D[u.Quantity]
    geometry: Geometry    

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

    @property
    def cell_centers_3d(self) -> "CellCenters3D":
        x1c_3d, x2c_3d, x3c_3d = np.meshgrid(self.x1c, self.x2c, self.x3c, indexing="ij")
        return CellCenters3D(
            x1c=x1c_3d,
            x2c=x2c_3d,
            x3c=x3c_3d,
            geometry=self.geometry,
        )

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

    def nH2(self, config:dict):
        MUSTAR = config["stars"]["mu_star"]
        numberrho_H2 = self.rho/(MUSTAR*uc.m_p.to(u.g))
        if numberrho_H2.unit!=1/u.cm**3:
            raise ValueError(f"number density unit should be 1/cm^3, not {numberrho_H2.unit}.")
        return numberrho_H2

#TODO: Opacity to be improved
@dataclass(kw_only=True, slots=True)
class Opacity:
    mix:dict
    rho:float|None=None
    value:str=field(init=False) # not always a string, can be instance of do.diel_const

    def __post_init__(self):
        if self.mix["mode"]=="file":
            if "reference" in self.mix:
                print(f"Please cite {self.mix['reference']} when using these optical constants.")
            self.value = do.diel_from_lnk_file(
                self.mix["file"], 
                headerlines=self.mix["headerlines"], 
            )
            if not is_set(self.rho):
                raise ValueError(
                    f"Internal density of the mix has to be defined. Please provide 'rho' in dust.opacity."
                )
            self.value.rho = self.rho
            set_lambda = set(np.sign(np.diff(self.value._l)))
            if len(set_lambda)!=1:
                raise ValueError(
                    f"Optical constants should be ordered by monotonically increasing lambda"
                )
            if list(set_lambda)[0]==-1:
                print("WARNING: lambda is monotically decreasing. Reversing optical constants arrays.")
                value_dict = vars(self.value)
                for key in ("_l","_n","_k","_ll","_ln","_lk"):
                    value_dict[key] = value_dict[key][::-1]
            if extrapolate:=self.mix["extrapolate_lambda_micron"]:
                mandatory_extrapolate_keys = {"min","max","N"}
                if extrapolate["mode"] in ("up","down") and set(extrapolate.keys())-{"mode"}!=mandatory_extrapolate_keys:
                    raise ValueError(
                        f"{(set(extrapolate.keys())-{'mode'}) ^ mandatory_extrapolate_keys} should be specified in 'extrapolate_lambda_micron'."
                    )
                lmin = (extrapolate["min"]*u.micron).to(u.cm).value
                lmax = (extrapolate["max"]*u.micron).to(u.cm).value
                if extrapolate["mode"]=="up":
                    self.value.extrapolate_constants_up(
                        lmin=lmin, 
                        lmax=lmax, 
                        n=extrapolate["N"], 
                        kind="second",
                    )
                elif extrapolate["mode"]=="down":
                    self.value.extrapolate_constants_down(
                        lmin=lmin, 
                        lmax=lmax, 
                        n=extrapolate["N"], 
                        kind="second",
                    )
                else:
                    raise ValueError(
                        f"Unknown extrapolation mode in 'extrapolate_lambda_micron': {extrapolate['mode']}. Should be 'up' or 'down'."
                    )

        elif self.mix["mode"] in ("birnstiel2018","ricci2010"):
            self.value = self.mix["mode"]
            if is_set(self.rho):
                print("WARNING: unused 'rho' when using dsharp_opac mix.")
        else:
            raise ValueError(
                f"Unknown mode for dust opacity mix: {self.mix['mode']}. Should be 'file' or 'dsharp_opac'."
            )        

@dataclass(kw_only=True, slots=True)
class Model:
    grid:Grid|None=None
    dust:Dust|None=None
    gas:Gas|None=None
    unit_length_au:float
    unit_mass_msun:float
    component:str

    @property
    def dimension(self) -> int:
        return int(3-self.grid.shape.count(1))

    @property
    def axes_from_geometry(self):
        return axes_from_geometry(self.grid.geometry)

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

    def temperature(self, config:dict) -> "Array":
        if "gas" not in self.component:
            raise ValueError(
                "temperature method is defined if there is a gas component."
            )
        config_gas = config["gas"]
        if temperature_dict:=config_gas["temperature"]:
            mode = temperature_dict["mode"]
        else:
            raise ValueError(
                "Undefined treatment of temperature. Either use 'tgas_eq_tdust = 1' in radmc3d.inp, or provide 'temperature' in the gas section."
            )
        match mode:
            case "array":
                if "processing" not in config["gas"]:
                    raise ValueError(
                        f"'processing' has to be provided in the gas section in order to retrieve the temperature array from prodimo."
                    )
                processing_dict = np.atleast_1d(config["gas"]["processing"]).tolist()
                processing_category = [d.get("mode") for d in processing_dict]
                if len(processing_category)!=1 and "prodimo" not in processing_category:
                    raise NotImplementedError(
                        f"For now, we can only retrieve temperature array with prodimo mode. No other mode is implemented yet."
                    )
                if "prodimo_dir" not in processing_dict[0]:
                    raise ValueError(
                        f"'prodimo_dir' should be provided, when using mode='prodimo'."
                    ) 
                prodimo_directory = processing_dict[0]["prodimo_dir"]
                pmodel = pread.read_prodimo(prodimo_directory)

                x_2d, z_2d = ((pmodel.x*u.au).to(u.cm), (pmodel.z*u.au).to(u.cm))

                x_3d = np.atleast_3d(x_2d)
                z_3d = np.atleast_3d(z_2d)
                pcyl_3d = np.zeros_like(x_3d.value)*u.radian

                temperature = Array(
                    cells=CellCenters3D(
                        x1c=x_3d,
                        x2c=z_3d,
                        x3c=pcyl_3d,
                        geometry=Geometry("polar"),
                    ),
                    data=np.atleast_3d(pmodel.tg)*u.K,
                    config=config,
                )
                temperature = temperature.from_prodimo_to(output_cells=self.grid.cell_centers_3d)
                if temperature.data.shape!=self.gas.rho.shape:
                    raise ValueError(
                        f"temperature and gas density do not have the same shape: {temperature.data.shape}!={self.gas.rho.shape}. "\
                        "Try to perform some array processing."
                    )
                if temperature.data.unit!=u.K:
                    raise ValueError(f"temperature unit should be K, not {temperature.data.unit}.")

            case _ as unreachable:
                assert_never(unreachable)
        return temperature

    def number_density(self, config:dict) -> "Array":
        if "gas" not in self.component:
            raise ValueError(
                "number_density method is defined if there is a gas component."
            )
        config_gas = config["gas"]
        species = config_gas["species"]
        abundance = config_gas["abundance"]
        mode = abundance["mode"]
        value = abundance["value"]

        match mode, value:
            case "unset", "unset":
                if "processing" not in config["gas"]:
                    raise ValueError(
                        f"'processing' has to be provided in the gas section in order to retrieve the number density array from prodimo."
                    )
                processing_dict = np.atleast_1d(config["gas"]["processing"]).tolist()
                processing_category = [d.get("mode") for d in processing_dict]
                if len(processing_category)!=1 and "prodimo" not in processing_category:
                    raise NotImplementedError(
                        f"For now, we can only retrieve number density array with prodimo mode. No other mode is implemented yet."
                    )
                if "prodimo_dir" not in processing_dict[0]:
                    raise ValueError(
                        f"'prodimo_dir' should be provided, when using mode='prodimo'."
                    ) 
                prodimo_directory = processing_dict[0]["prodimo_dir"]
                pmodel = pread.read_prodimo(prodimo_directory)

                x_2d, z_2d = ((pmodel.x*u.au).to(u.cm), (pmodel.z*u.au).to(u.cm))

                x_3d = np.atleast_3d(x_2d)
                z_3d = np.atleast_3d(z_2d)
                pcyl_3d = np.zeros_like(x_3d.value)*u.radian

                distribution = Array(
                    cells=CellCenters3D(
                        x1c=x_3d,
                        x2c=z_3d,
                        x3c=pcyl_3d,
                        geometry=Geometry("polar"),
                    ),
                    data=np.atleast_3d(pmodel.nmol[:,:,pmodel.spnames[species.upper()]])*(1/u.cm**3),
                    config=config,
                )
                distribution = distribution.from_prodimo_to(output_cells=self.grid.cell_centers_3d)
                nH2 = self.gas.nH2(config=config)
                if distribution.data.shape!=nH2.shape:
                    raise ValueError(
                        f"number density of species and nH2 do not have the same shape: {distribution.data.shape}!={nH2.shape}. "\
                        "Try to perform some array processing."
                    )
                if distribution.data.unit!=1/u.cm**3:
                    raise ValueError(f"number density unit should be 1/cm^3, not {distribution.data.unit}.")

            case "array", "unset":
                if "processing" not in config["gas"]:
                    raise ValueError(
                        f"'processing' has to be provided in the gas section in order to retrieve the number density array from prodimo."
                    )
                processing_dict = np.atleast_1d(config["gas"]["processing"]).tolist()
                processing_category = [d.get("mode") for d in processing_dict]
                if len(processing_category)!=1 and "prodimo" not in processing_category:
                    raise NotImplementedError(
                        f"For now, we can only retrieve number density array with prodimo mode. No other mode is implemented yet."
                    )
                if "prodimo_dir" not in processing_dict[0]:
                    raise ValueError(
                        f"'prodimo_dir' should be provided, when using mode='prodimo'."
                    ) 
                prodimo_directory = processing_dict[0]["prodimo_dir"]
                pmodel = pread.read_prodimo(prodimo_directory)

                x_2d, z_2d = ((pmodel.x*u.au).to(u.cm), (pmodel.z*u.au).to(u.cm))

                x_3d = np.atleast_3d(x_2d)
                z_3d = np.atleast_3d(z_2d)
                pcyl_3d = np.zeros_like(x_3d.value)*u.radian

                abundance = Array(
                    cells=CellCenters3D(
                        x1c=x_3d,
                        x2c=z_3d,
                        x3c=pcyl_3d,
                        geometry=Geometry("polar"),
                    ),
                    data=np.atleast_3d(pmodel.getAbun(species.upper()))*u.dimensionless_unscaled,
                    config=config,
                )
                abundance = abundance.from_prodimo_to(output_cells=self.grid.cell_centers_3d)
                nH2 = self.gas.nH2(config=config)
                if abundance.data.shape!=nH2.shape:
                    raise ValueError(
                        f"abundance and nH2 do not have the same shape: {abundance.data.shape}!={nH2.shape}. "\
                        "Try to perform some array processing."
                    )
                distribution = Array(
                    cells=abundance.cells,
                    data=abundance.data*nH2,
                    config=config,
                )
                if distribution.data.unit!=1/u.cm**3:
                    raise ValueError(f"number density unit should be 1/cm^3, not {distribution.data.unit}.")

            case "constant", float:
                abundance = value
                nH2 = self.gas.nH2(config=config)

                distribution = Array(
                    cells=CellCenters3D(
                        x1c=self.grid.cell_centers_3d.x1c,
                        x2c=self.grid.cell_centers_3d.x2c,
                        x3c=self.grid.cell_centers_3d.x3c,
                        geometry=self.grid.cell_centers_3d.geometry,
                    ),
                    data=abundance*nH2,
                    config=config,
                )
            case _ as unreachable:
                assert_never(unreachable)
        return distribution

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
                    opacity=opacity.value,
                    smoothing=smoothing,
                    config=config,
                )
            self._write_metadata(directory=directory)
        if "gas" in self.component:
            self._write_numberdens_inp(
                directory=directory, 
                config=config
            )
            self._write_gas_velocity_inp(directory=directory)
            if config["gas"]["temperature"]:
                self._write_gastemperature_inp(
                    directory=directory, 
                    config=config
                )

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
        ### Nspec = 100
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
            ### a_opac = np.geomspace(dustsize.min(), dustsize.max(), Nspec)
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
                        f"{opac_dict['k_abs'][ia, ilam]:.6e} "\
                        f"{opac_dict['k_sca'][ia, ilam]:.6e} "\
                        f"{opac_dict['g'][ia, ilam]:.6e}\n"
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
        species = config["gas"]["species"]
        numberrho = self.number_density(config=config).data

        filename = f"numberdens_{species}.binp"
        path = os.path.join(directory, filename)

        print(f"INFO: Writing {path}.....", end="")
        with open(path, "wb") as f:
            header = np.array([1, 8, numberrho.flatten().shape[0]], dtype=int)
            header.tofile(f)
            numberrho.ravel(order="F").value.tofile(f)
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
        print("done.")

    def _write_gastemperature_inp(self, *, directory:str, config:dict):
        """
        Function writes the 'gas_temperature.inp' input file.

        Parameters
        ----------
        directory : str
            Data directory in which the files are written.
        """
        temperature = self.temperature(config=config).data

        filename = f"gas_temperature.binp"
        path = os.path.join(directory, filename)

        print(f"INFO: Writing {path}.....", end="")
        with open(path, "wb") as f:
            header = np.array([1, 8, temperature.flatten().shape[0]], dtype=int)
            header.tofile(f)
            temperature.ravel(order="F").value.tofile(f)
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
        if self.grid.geometry!="spherical":
            raise NotImplementedError(f"Grid geometry={self.grid.geometry} should be 'spherical'.")
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
    )
    return model
