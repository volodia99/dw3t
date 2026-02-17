# dw3t
[![PyPI](https://img.shields.io/pypi/v/dw3t.svg?logo=pypi&logoColor=white&label=PyPI)](https://pypi.org/project/dw3t/)
[![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://github.com/astral-sh/uv)

***Interface between spherical datasets from simulation outputs (e.g., readable with [nonos](https://github.com/la-niche/nonos)) to radmc3d.***

## Development status

*Word of caution*: `dw3t` has still to be tested, in particular to be sure that everything works smoothly in a radmc3d computation. The documentation is work in progress.

***TODO***
- check optools ? astroquery ? instead of dsharp_opac
- workspace (1 repo, 2 distributions), for the interface radmc3d/prodimo or even a package (see pre-commit, idefix-cli)

## Installation

We recommend to install `dw3t` using the package and project manager `uv`. See the [documentation](https://docs.astral.sh/uv/getting-started/installation/#standalone-installer) to install `uv` on your system. Run the following:

```shell
uv tool install dw3t
```

## Use the interface

You can use the interface inside the project's virtual environment using a parameter .toml file: 

```shell
dw3t dw3t.toml
```

See the [TOML documentation](https://toml.io/en) to know more about this config file format.

## Configuration file

### Example

You can find an example for the parameter file in `dw3t/dw3t.toml`.

### 1. Section `[simulation]`

***Mandatory parameters:***
- `unit_length_au` : code unit of length \[au\] (`float`)
- `unit_mass_msun` : code unit of mass \[solMass\] (`float`)
- `component` : which component is included (`"dust"` and/or `"gas"`) (`str|list[str]`)

***Optional parameters:***
- `output_dir` : directory where the radmc3d files are computed (`str`). Default: `f"{input_dir}/radmc3d"`.
- `simulation_files_only` : computes only the radmc3d files related to the simulation, not the generic mandatory radmc3d files (`bool`). Default: `False`.
- `processing` : post-processing of the simulation dataset. The operations are chained from the first to the last dict-like element in the provided list. The final post-processed model has to be 3D (`list[dict]`). Default: no processing. There are two possible ways to use the `processing` parameter: 

### (i) with the templates

***Example:*** 
```toml
processing = [
    {mode="builtin", input_number=0, input_dir="tests/data/idefix_1_dust_fluid", internal_rho=2.0},
    {mode="phi_expansion", nphi=128},
]
```

Each element of the list has the form `{mode="choose_our_mode", **kwargs}`. For now, implemented modes are:
- "identity": returns a copy of the model (mainly for testing purposes)
- "builtin": use [nonos](https://github.com/la-niche/nonos) to read the grid and the fields necessary to compute a radmc3d model (tested only with idefix simulations for now). This mode comes necessary with `input_number` (`int`) and `input_dir` (`str`), which are respectively the number and the directory of the simulated output. Moreover, for idefix simulations you need to add the `internal_rho` (`float`) parameter (\[g/cm3\]) if `"dust"` in included in `component`.
- "phi_expansion": extend azimuthally a 2D spherical ($r$, $\theta$) simulation, using `nphi` cells in the $\phi$ direction, and returns a 3D spherical model.

See in `src/dw3t/template/model` for more info.

### (ii) with a user-defined python file

***Example:*** 
```toml
processing = {
    mode = "userdef",
    file = "src/dw3t/template/model/userdef.py",
    input_number = 0,
    input_dir = "tests/data/idefix_1_dust_fluid",
    internal_rho = 2.0,
}
```
When `mode = "userdef"`, only the argument `file` is necessary, which corresponds to the absolute path of the user-defined python processing file, which has to contain the `processing` function, with signature `processing(*, model:"Model", kwargs:dict) -> "Model"`, and returns a model object. This mode is for now incompatible with all other modes, so you need to post-process your data yourself in the given `file`.

See in `src/dw3t/template/model/userdef.py` for more info.

### 2. Section `[dust]`

This section for now contains only the information about the `opacity`, especially what dust mix is used, with its corresponding internal density if needed. We use the [dsharp_opac](https://github.com/birnstiel/dsharp_opac/) library to compute the optical constants associated to the mix. Two solutions are possible:

### (i) Using the available mixes in dsharp_opac

```toml
[dust.opacity]
mix = "birnstiel2018"
```
or
```toml
[dust.opacity]
mix = "ricci2010"
```

See the corresponding papers for more info.


### (ii) Using a .lnk file

Including the opacity constants. In that case you need to provide the internal density (`rho` \[g/cm3\]) of the mix.
```toml
[dust.opacity]
mix = "path_to_lnk_file"
rho = 2.0
```

### 3. Section `[gas]`

***Mandatory parameters:***

The two main parameters here are `species` and `number_density`. 
- `species`: name of the species (`str`). Used to construct the lines.inp, molecule_\*.dat and numberdens_\*.inp radmc3d files.
- `number_density`: computes the number density of the corresponding `species` (`dict`). Needs the arguments `abundance` (`dict`) (when defining a molecular abundance), and/or `processing` (`list[dict]`) depending on the need. The `abundance` dictionary is defined with the keys `mode` (`str`), and `value` (`float`).

***CASE 1***

When `mode="constant"`, the corresponding `value` must be defined:
```toml
[gas]
species = "co"

[gas.number_density]
abundance = {
    mode = "constant",
    value = 1e-4,
}
```
In that example, the gas number density resulting from the simulated output is multiplied by a constant abundance 1e-4, corresponding to a prescribed abundance for CO.

***CASE 2***

By default, an `abundance` that is not provided is unset. `processing` has then to be given to indicate how to retrieve the number density array associated to `species`, instead of computing it from the simulated gas density times `abundance`.
```toml
[gas]
species = "co"

[gas.number_density]
processing = [
    {mode="prodimo", input_dir="path_to_prodimo_model_directory"},
    {mode="midplane_symmetry"},
    {mode="phi_expansion", nphi=128},
]
```
In that example, the CO number density is retrieved from a 2D (prodimo)[https://prodimo.iwf.oeaw.ac.at] model. Th 2D array is then extended vertically by mirror symmetry and in 3D with an azimuthal expansion.

***CASE 3***

When `mode="array"`, `value` is `None` and the `processing` argument must be specified to indicate how to retrieve the abundance array associated to `species`.
```toml
[gas]
species = "co"

[gas.number_density]
abundance = {
    mode = "array",
}
```
***Remark:*** In CASE 2, we compute the CO number density array, whereas in CASE 3 we compute the CO abundance.

### 4. radmc3d-specific sections

There are 3 more sections, respectively linked to 3 radmc3d-related files : radmc3d.inp, stars.inp and wavelength_micron.inp.

#### (i) `[stars]` (mandatory):

Used in stars.inp. 

***Mandatory parameters:***
- `R_star`: radius of the star \[solRadius\] (`float`) 
- `T_star`: effective temperature of the star \[K\] (`float`). For now, we assume one unique star in the center of the grid. 

***Optional parameters:***
- `M_star`: mass of the star \[solMass\] (`float`) default: `unit_mass_msun`.
- `mu_star`: mean molecular weight of the star (`float`) default: `1.37`.

Example:
```toml
[stars]
R_star = 2.0
T_star = 7_000
nphot = 1_000
nphot_scat = 1_000
```

#### (ii) `[radmc3d]` (optional):

Used in radmc3d.inp, with the same parameter names as the ones defined in radmc3d. Example:
```toml
[radmc3d]
istar_sphere = 1
tgas_eq_tdust = 1
nphot = 1_000
nphot_scat = 1_000
```

***Careful: by default, we construct an empty radmc3d.inp file.***

#### (iii) `[wavelength_micron]` (optional):

Used in stars.inp, wavelength_micron.inp. Corresponds to the wavelength coverage used to define the stellar spectrum, from `min` (`float`) to `max` (`float`) with `N` (`int`) wavelength points. Note that for simplicity we use the same definition for the wavelength coverage to construct the dustkapscatmat_\*.inp file. By default:
```toml
[wavelength_micron]
min = 1
max = 10_000
N = 200
```
