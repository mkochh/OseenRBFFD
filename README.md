# OseenRBFFD
## Installation
This project requires the CGAL library ([https://www.cgal.org/](https://www.cgal.org/)). 
The Code should work with CGAL 6.x and 5.x (although tested only on 6.0.1 and 5.5.1), 
it might also work with older versions, but this was not tested.

To install CGAL and other depencies on Ubuntu:
```
sudo apt-get install git g++ python3 cmake libhdf5-serial-dev doxygen graphviz libcgal-dev libopenblas-dev
```

or on OpenSuse:
```
sudo zypper install git gcc python3 cmake hdf5-devel doxygen graphviz cgal-devel openblas-devel
```

Next clone this repository:
```
git clone https://github.com/mkochh/OseenRBFFD
```

## Compilation
Execute the follwing from your terminal to compile the code:
```
cd OseenRBFFD
cmake .
make
```

Wait for compilation to finish. This might take some time when compiling the code for the first time.
You can use `make -j n` to allow make to use parallelization and speed up the process, where `n` is the number of threads to parallize on.
Beware that this can be quite ressource intensive, if `n` is close to your maximal number of threads.

If you get compile errors regarding HDF5, please check [https://e6.ijs.si/medusa/wiki/index.php/How_to_build#HDF5](https://e6.ijs.si/medusa/wiki/index.php/How_to_build#HDF5) for a guide on how to fix this.

## Check that everything works:
Execute the following form your terminal to run the executable `oseenHLU`:
```
cd Tests
./oseenHLU
```
You should see some output informing you about the progress of the program.
Once it is finished you should find a .csv file in OseenRBFFD/Tests/Daten.
