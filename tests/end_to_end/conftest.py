import pytest
from pathlib import Path
import subprocess
import sys
from typing import Callable


@pytest.fixture
def output_dir() -> str:
    return "test_output"


@pytest.fixture
def run() -> Callable:
    def func(config_file: Path):
        command_list = ["./earlGreyParTEA", "--threads", "4", "--config", config_file]
        process = subprocess.Popen(command_list)
        process.wait()
    return func
