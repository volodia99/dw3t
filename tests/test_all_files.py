import os
from pathlib import Path
import tarfile
import filecmp
import shutil
import pytest

from dw3t.main import main

def partial_extraction(*, target_tarfile:str, target_subdirectory:str, path:str):
    with tarfile.open(target_tarfile) as tar:
        subdir_and_files = [
            tarinfo for tarinfo in tar.getmembers()
            if tarinfo.name.startswith(target_subdirectory)
        ]
        tar.extractall(path=path, members=subdir_and_files, filter="data")

class TestFileWrite:
    def test_all_files(self, test_data_dir, data_dir):
        if os.path.isdir(test_data_dir / "idefix_1_dust_fluid"):
            shutil.rmtree(test_data_dir / "idefix_1_dust_fluid")
        partial_extraction(
            target_tarfile=data_dir / "idefix_1_dust_fluid.tar.gz", 
            target_subdirectory="idefix_1_dust_fluid/radmc3d_ref/",
            path=test_data_dir,
        )
        partial_extraction(
            target_tarfile=data_dir / "idefix_1_dust_fluid.tar.gz", 
            target_subdirectory="idefix_1_dust_fluid/idefix_ref/",
            path=test_data_dir,
        )
        data_dir_radmc3d_ref = test_data_dir / "idefix_1_dust_fluid" / "radmc3d_ref"
        data_dir_radmc3d = test_data_dir / "idefix_1_dust_fluid" / "radmc3d"

        main([str(test_data_dir / "dw3t.toml")])

        assert filecmp.cmp(data_dir_radmc3d_ref / "amr_grid.inp", data_dir_radmc3d / "amr_grid.inp", shallow=False)
        assert filecmp.cmp(data_dir_radmc3d_ref / "dust_density.binp", data_dir_radmc3d / "dust_density.binp", shallow=False)
        for kk in range(1):
            assert filecmp.cmp(data_dir_radmc3d_ref / f"dustkapscatmat_{kk}.inp", data_dir_radmc3d / f"dustkapscatmat_{kk}.inp", shallow=False)
        assert filecmp.cmp(data_dir_radmc3d_ref / "dustopac.inp", data_dir_radmc3d / "dustopac.inp", shallow=False)
        assert filecmp.cmp(data_dir_radmc3d_ref / "gas_velocity.binp", data_dir_radmc3d / "gas_velocity.binp", shallow=False)
        assert filecmp.cmp(data_dir_radmc3d_ref / "lines.inp", data_dir_radmc3d / "lines.inp", shallow=False)
        assert filecmp.cmp(data_dir_radmc3d_ref / "metadata.dw3t.npz", data_dir_radmc3d / "metadata.dw3t.npz", shallow=False)
        assert filecmp.cmp(data_dir_radmc3d_ref / "molecule_co.inp", data_dir_radmc3d / "molecule_co.inp", shallow=False)
        assert filecmp.cmp(data_dir_radmc3d_ref / "numberdens_co.binp", data_dir_radmc3d / "numberdens_co.binp", shallow=False)
        assert filecmp.cmp(data_dir_radmc3d_ref / "radmc3d.inp", data_dir_radmc3d / "radmc3d.inp", shallow=False)
        assert filecmp.cmp(data_dir_radmc3d_ref / "stars.inp", data_dir_radmc3d / "stars.inp", shallow=False)
        assert filecmp.cmp(data_dir_radmc3d_ref / "wavelength_micron.inp", data_dir_radmc3d / "wavelength_micron.inp", shallow=False)

        shutil.rmtree(test_data_dir / "idefix_1_dust_fluid")        