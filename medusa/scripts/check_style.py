#/usr/bin/python3
import os
import sys
import multiprocessing as mp
import subprocess as sp

cwd = os.getcwd()
folder = cwd.split('/')[-1]
assert folder == 'medusa', "Please run from the root medusa/ folder."

STYLEFILTERS = ','.join([
    "-legal",
    "-build/include",
    "-runtime/reference",
    "-runtime/threadsafe_fn",
    "-runtime/explicit",
    "-readability/streams",
    "-whitespace/empty_loop_body",
    "-build/c++11",
    "-runtime/int",
])

FOLDERS = ['include/medusa', 'src', 'test', 'examples']
ENDINGS = frozenset(['hpp', 'cpp'])
EXCLUDE_DIRS = frozenset(['build'])

to_process = []
for folder in FOLDERS:
    for (path, dirs, files) in os.walk(folder):
        i = 0
        while i < len(dirs):
            if dirs[i] in EXCLUDE_DIRS:
                del dirs[i]
            else:
                i += 1

        for file in files:
            end = file.split('.')[-1]
            if end in ENDINGS:
                to_process.append(os.path.join(path, file))

def check_style(filename):
    try:
        cmd = [sys.executable, 'scripts/cpplint.py',
                               '--filter='+STYLEFILTERS,
                               '--linelength=100',
                               '--root=include',
                               '--extensions=hpp,cpp',
                               filename]
        p = sp.Popen(cmd, stderr=sp.PIPE, stdout=sp.PIPE)
        (stdout, stderr) = p.communicate()
    except:
        p.kill()
        (stdout, stderr) = p.communicate()

    if p.returncode != 0:
        lines = stderr.decode().split('\n')
        print('\n'.join(lines[:-3]))
    return p.returncode

def run_stylechek():
    cc = mp.cpu_count()
    with mp.Pool(cc) as pool:
        errors = pool.map(check_style, to_process)
    return sum(errors)

if __name__ == '__main__':
    run_stylechek()
