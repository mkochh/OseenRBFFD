#!/bin/bash

# this script checks system' configuration

BW='\x1b[37;1m'  # bold white
BR='\x1b[31;1m'  # bold red
BG='\x1b[32;1m'  # bold red
NC='\x1b[37;0m'  # no color

# packages and commands -- default = ubuntu
HDF5="libhdf5-dev"  # libhdf5-cpp-8

CHECK_PACKAGE="dpkg -s"
INSTALL_PACKAGE="sudo apt-get install"

if command -v pacman > /dev/null 2>&1; then
    echo "Detected Arch-like disto."
    CHECK_PACKAGE="pacman -Q"
    INSTALL_PACKAGE="sudo pacman -S"
    HDF5="hdf5"
    IHDF5=""
fi

TO_INSTALL=""
FAILED=0

function error { >&2 echo -e "${BR}$1${NC}"; }
function showinfo { echo -e "${BW}$1${NC}"; }
function ok { echo -e "${BG}$1${NC}"; }

function check_installed {
    echo "    Checking $1 ..."
    $1 --version > /dev/null
    if [ $? -ne 0 ]; then
        error "You dont have $1 installed."
        TO_INSTALL+=" $1"
        FAILED=1
    fi
}

function check_package {
    for pac in "$@"; do
        echo "    Checking package $pac ..."
        $CHECK_PACKAGE $pac > /dev/null
        if [ $? -ne 0 ]; then
            error "You dont have $pac installed."
            echo -e "Please run\n    $INSTALL_PACKAGE $pac\nto install."
            TO_INSTALL+=" $pac"
            FAILED=1
        fi
    done
}

cwd=${PWD##*/}
if [ "$cwd" != "medusa" ]; then
    echo -e "${BR}Run from medusa/ directory in this project!${NC}"
    exit 1
fi

showinfo "Checking build tools ..."
check_installed cmake
check_installed make

showinfo "Checking package dependencies ..."
check_package $HDF5
if [ $FAILED -eq 1 ]; then
    error "Missing packages!"
    echo -e "Run"
    echo "    $INSTALL_PACKAGE $TO_INSTALL"
    exit 1
fi
