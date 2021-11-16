import os
import shutil
import tempfile

import pytest


@pytest.fixture(scope="class")
def cleandir():
    old_cwd = os.getcwd()
    newpath = tempfile.mktemp(dir=old_cwd)
    os.chdir(newpath)
    newfile = tempfile.mkstemp()
    yield newpath, newfile
    os.chdir(old_cwd)
    shutil.rmtree(newpath)
