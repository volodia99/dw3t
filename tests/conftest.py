from pathlib import Path

import pytest

@pytest.fixture()
def data_dir():
    return Path(Path(__file__).parent).parent / "datasets"

@pytest.fixture()
def test_data_dir():
    return Path(__file__).parent
