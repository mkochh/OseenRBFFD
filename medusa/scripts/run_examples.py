import multiprocessing as mp
import subprocess as sp
import os

def run_examples():
    cwd = os.getcwd()
    folder = cwd.split('/')[-1]
    assert folder == 'medusa', "Please run from the root medusa/ folder."

    ret = 0
    try:
        if os.path.exists('build') and not os.path.isdir('build'):
            print('build/ exists, but is not a directory. Please fix.')
            ret = 1
            return
        if not os.path.exists('build'):
            os.mkdir('build')
        os.chdir('build')
        ret = sp.call('cmake ..', shell=True)
        if ret != 0: return
        cc = mp.cpu_count()
        ret = sp.call('make -j {} examples_run'.format(cc), shell=True)
        if ret != 0: return
    finally:
        os.chdir(cwd)
    return ret

if __name__ == '__main__':
    run_examples()
