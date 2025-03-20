import os
import subprocess as sp

def configure():
    cwd = os.getcwd()
    folder = cwd.split('/')[-1]
    assert folder == 'medusa', "Please run from the root medusa/ folder."

    return sp.call('./scripts/configure.sh', shell=True)
