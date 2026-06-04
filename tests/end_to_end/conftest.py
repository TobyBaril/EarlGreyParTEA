import pytest
from pathlib import Path
import subprocess
from typing import Callable


@pytest.fixture
def output_dir() -> str:
    return "test_output"


@pytest.fixture
def run() -> Callable:
    def func(config_file: Path):
        command_list = ["./earlGreyParTEA", "--threads", "4", "--config", str(config_file)]
        process = subprocess.Popen(command_list)
        process.wait()
        if process.returncode != 0:
            raise RuntimeError(f"Command {' '.join(command_list)} exited with code {process.returncode}")
    return func
