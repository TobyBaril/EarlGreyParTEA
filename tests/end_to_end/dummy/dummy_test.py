import pytest
from pathlib import Path
from typing import Callable

import yaml


class TestDummy:

    def test_dummy_fasta_3_times(self, run: Callable):
        config_file = Path(__file__).parent / "dummy_config.yaml"
        run(config_file)
