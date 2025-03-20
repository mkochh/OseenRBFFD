# Compiling and running examples

Our standard example name will be `poisson_dirichlet_1D`, but the instructions naturally apply to all other
examples. Once you have built medusa, you can build and run the examples from your build
directory simply as
```
make poisson_dirichlet_1d_run
```
This will run the example from its own directory, e.g. `examples/poisson_equation`.

# Standalone compilation
Alternatively, you may wish to compile and run the examples individually. First, you have to build
the `medusa_standalone` target in your medusa build directory, which creates a static library
`medusa_standalone`.

The examples can then be compiled manually using e.g. `g++`, provided that
appropriate include and library paths are set and the `medusa_standalone` library is specified
for linking.

Example command to be run from `poisson_equation` directory:
```
g++ -o poisson_dirichlet_1D poisson_dirichlet_1D.cpp -I ../../include -Wall -O3 -L ../../bin -lmedusa_standalone
```

The compiled example can be run simply as
```
./poisson_dirichlet_1D
```

More detailed instuctions and the description of the philosphy of the examples is available on
[our wiki](http://e6.ijs.si/medusa/wiki/index.php/Philosophy_of_examples_and_how_to_run_them).
