DEFAULT_LAYER = {
    "simulation": {
        "output_dir": "unset",
        "simulation_files_only": False,
        "processing": {
            "mode": "identity",
        },
    },
    "dust": {
        "opacity": {
            "mix": {
                "extrapolate_lambda_micron": {},
            },
            "rho": "unset",
            "smoothing": False,
        },
    },
    "radmc3d": {},
    "stars": {
        "M_star": "unset",
        "mu_star": 1.37,
    },
    "wavelength_micron": {
        "min": 1,
        "max": 10_000,
        "N": 200,
    },
}

MANDATORY_SET = {
    "unit_length_au",
    "unit_mass_msun",
    "component",
    "R_star",
    "T_star",
}
