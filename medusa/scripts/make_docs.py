import subprocess as sp
import os

def make_docs():
    cwd = os.getcwd()
    folder = cwd.split('/')[-1]
    assert folder == 'medusa', "Please run from the root medusa/ folder."
    try:
        p = sp.Popen(['doxygen', 'docs/Doxyfile'], stdout=sp.PIPE, stderr=sp.PIPE)
        (stdout, stderr) = p.communicate()
    except:
        p.kill()
        (stdout, stderr) = p.communicate()

    if stderr:
        print(stderr.decode(), end='')
        return 1
    return 0

if __name__ == '__main__':
    make_docs()
