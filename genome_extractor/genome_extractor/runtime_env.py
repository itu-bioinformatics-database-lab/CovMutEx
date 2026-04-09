import os
import tempfile
from pathlib import Path


def configure_runtime_tempdir():
    project_dir = Path(__file__).resolve().parent.parent
    temp_dir = project_dir / ".tmp"
    temp_dir.mkdir(parents=True, exist_ok=True)

    temp_dir_str = str(temp_dir)
    for env_name in ("TMPDIR", "TEMP", "TMP", "TEST_TMPDIR"):
        os.environ[env_name] = temp_dir_str

    tempfile.tempdir = temp_dir_str
    return temp_dir_str
