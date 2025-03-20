# OseenRBFFD
## Installation
This project requires the CGAL library ([https://www.cgal.org/](https://www.cgal.org/)). 
Therefore, install CGAL (libcgal-dev or cgal-devel or similar depending on your OS with a packet manager).
The Code should work with CGAL 6.x and 5.x (although tested only on 6.0.1 and 5.5.1), 
it might also work with older versions, but this was not tested.

Next clone this repository:
```
git clone https://github.com/mkochh/OseenRBFFD
```

## Compilation
You will need cmake and gcc.
In the directory OseenRBFFD, run:
```
cmake .
make -j
```

Wait for compilation to finish. This might take some time when compiling the code for the first time.

## Check that everything works:
Go in to the directory Tests and run the executable `oseenHLU`:
```
cd Tests
./oseenHLU
```
You should see some output informing you about the progress of the program.
Once it is finished you should find a .csv file in OseenRBFFD/Tests/Daten.
