import os
import filecmp
import shutil
import pytest

from dw3t.main import main

class TestFileWrite:
    def test_all_files(self, test_data_dir):
        data_dir_ref = test_data_dir / "idefix_6_dust_fluids" / "radmc3d_ref"
        data_dir = test_data_dir / "idefix_6_dust_fluids" / "radmc3d"

        if os.path.isdir(data_dir):
            shutil.rmtree(data_dir)
        main([str(test_data_dir / "dw3t.ini")])

        assert filecmp.cmp(data_dir_ref / "amr_grid.inp", data_dir / "amr_grid.inp", shallow=False)
        assert filecmp.cmp(data_dir_ref / "dust_density.binp", data_dir / "dust_density.binp", shallow=False)
        for kk in range(6):
            assert filecmp.cmp(data_dir_ref / f"dustkapscatmat_{kk}.inp", data_dir / f"dustkapscatmat_{kk}.inp", shallow=False)
        assert filecmp.cmp(data_dir_ref / "dustopac.inp", data_dir / "dustopac.inp", shallow=False)
        assert filecmp.cmp(data_dir_ref / "gas_velocity.binp", data_dir / "gas_velocity.binp", shallow=False)
        assert filecmp.cmp(data_dir_ref / "lines.inp", data_dir / "lines.inp", shallow=False)
        assert filecmp.cmp(data_dir_ref / "metadata.npz", data_dir / "metadata.npz", shallow=False)
        assert filecmp.cmp(data_dir_ref / "molecule_co.inp", data_dir / "molecule_co.inp", shallow=False)
        assert filecmp.cmp(data_dir_ref / "numberdens_co.binp", data_dir / "numberdens_co.binp", shallow=False)
        assert filecmp.cmp(data_dir_ref / "radmc3d.inp", data_dir / "radmc3d.inp", shallow=False)
        assert filecmp.cmp(data_dir_ref / "stars.inp", data_dir / "stars.inp", shallow=False)
        assert filecmp.cmp(data_dir_ref / "wavelength_micron.inp", data_dir / "wavelength_micron.inp", shallow=False)

        shutil.rmtree(test_data_dir / "idefix_6_dust_fluids" / "radmc3d")