#!/usr/bin/python3
import argparse
parser = argparse.ArgumentParser(description='Script to check you medusa library.')
parser.add_argument('-c', '--configure', action='store_true', help='check your system configuration')
parser.add_argument('-b', '--build', action='store_true', help='build the library')
parser.add_argument('-t', '--tests', action='store_true', help='run all tests')
parser.add_argument(      '--build-examples', action='store_true', help='build all examples')
parser.add_argument('-e', '--run-examples', action='store_true', help='run selected examples')
parser.add_argument('-s', '--style', action='store_true', help=(
    'style checks all sources and reports any errors'))
parser.add_argument('-d', '--docs', action='store_true', help=(
    'generate and check for documentation errors'))

args = parser.parse_args()
if not any(args.__dict__.values()):
    args.tests = True
    args.build_examples = True
    args.run_examples = True
    args.style = True
    args.docs = True

if args.tests: args.build = True
if args.run_examples: args.build_examples = True

import os
folder = os.getcwd().split('/')[-1]
assert folder == 'medusa', "Please run from the root medusa/ folder."

class CC:
    BOLD_WHITE ='\x1b[37;1m'  # bold white
    BOLD_RED ='\x1b[31;1m'  # bold red
    BOLD_GREEN = '\x1b[32;1m'  # bold green
    NO_COLOR = '\x1b[37;0m'  # no color

def info(text):
    print(CC.BOLD_WHITE+text+CC.NO_COLOR)

def error(text):
    print(CC.BOLD_RED+text+CC.NO_COLOR)

def good(text):
    print(CC.BOLD_GREEN+text+CC.NO_COLOR)


def main():
    if args.configure:
        info('Configuring...')
        import scripts.configure
        ret = scripts.configure.configure()
        if ret != 0:
            error("Your system is missing critical dependencies!")
            return False

    if args.build:
        info("Building...")
        import scripts.build
        ret = scripts.build.build()
        if ret != 0:
            error("There were errors during building!")
            return False

    if args.tests:
        info("Running tests...")
        import subprocess as sp
        ret = sp.call('./bin/medusa_tests')
        if ret != 0:
            error("There were failed tests.")
            return False

    if args.style:
        info("Checking style...")
        import scripts.check_style
        num_errors = scripts.check_style.run_stylechek()
        if num_errors != 0:
            error("There were style errors in {} files!".format(num_errors))
            return False

    if args.docs:
        info("Generating documentation...")
        import scripts.make_docs
        errors = scripts.make_docs.make_docs()
        if errors != 0:
            error("There were documentation errors!")
            return False

    if args.build_examples:
        info("Building examples...")
        import scripts.build_examples
        ret = scripts.build_examples.build_examples()
        if ret != 0:
            error("There were errors when building examples!")
            return False

    if args.run_examples:
        info("Running examples...")
        import scripts.run_examples
        ret = scripts.run_examples.run_examples()
        if ret != 0:
            error("There were errors when running examples!")
            return False

    return True

if __name__ == '__main__':
    import sys
    ret = main()
    if ret == True:
        good('Everything looks good.')
        sys.exit(0)
    sys.exit(1)
