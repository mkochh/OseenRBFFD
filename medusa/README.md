# Medusa

Coordinate Free Meshless Method implementation.
See http://e6.ijs.si/medusa for more details.

## Installation
To make this work from plain Ubuntu installation, run
```
sudo apt-get install git g++ python3 cmake libhdf5-dev doxygen graphviz
git clone https://gitlab.com/e62Lab/medusa.git --branch master --single-branch
cd medusa
./run_tests.py
```
which installs dependencies, clones the repository, goes into the root folder of
the repository and runs tests. This will build and run all tests. If this
works, you are ready to go! Otherwise install any missing packages and if it
still fails, raise an issue.

For building on other systems and troubleshooting refer to the
[Installation and building guide](http://e6.ijs.si/medusa/wiki/index.php/How_to_build).

For instructions on how to use this library in you project, see
[Including this library in your project](http://e6.ijs.si/medusa/wiki/index.php/Including_this_library_in_your_project).

E62Lab team
