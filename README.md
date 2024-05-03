# Eigenious: Linear Algebra Library

Eigenious is a project involving writing a C++ library for linear algebra and matrix computations, including solving eigenvalue problems. Within the library are two classes:
- `Matrix`: A class for storing and manipulating matrices
- `Vector`: A class for storing and manipulating vectors

Both of these classes are templated, meaning that they can be used to store any type of data. 
However, the library is primarily intended for use with floating point numbers, and is to be used only with int, float, and double types for numerical computations.

Diagonalization of the matrix is done using the Jacobi method, which is an iterative method for finding the eigenvalues and eigenvectors of a matrix involving repeated rotations of the matrix until it is diagonalized. The implementation of the Jacobi method is found in the `Matrix` class.

All of the functions in the library are tested using cassert, and the tests can be found in the `test` folder. The library is also tested on two Quantum Mechanics problem sets: an extended Huckle program and CNDO/2 program. These programs can be found in the `Tests` folder
within 'ProblemSets'.

## Creating the Library and Running the Tests
There is a makefile included in the project, which can be used to compile the library and run the tests. To compile the library, run `make all` in the terminal. 

The library `eigenious.a` will be created in the `Lib` folder, and the tests will be compiled into the `bin` folder. 

To run the tests, enter the `Bin` folder and type `./` followed by one of the test names.

The following tests do not require any input files:
- `test_matrix`: Tests the `Matrix` class
- `test_vector`: Tests the `Vector` class
- `test_matrix_ops`: Tests the matrix operations in the `Matrix` class
- `test_vector_ops`: Tests the vector operations in the `Vector` class
- `ps3_energy_diff`: Tests the extended Huckle program on a specific example
- `ps4_oxygen`: Tests the CNDO/2 program on oxygen

The following tests require input files found in the `Molecules` directory:
- `ps3_test`: Tests the extended Huckle program on the input molecule
- `ps4_test`: Tests the CNDO/2 program on the input molecule

## Cleaning the project
To clean the project, run `make clean` in the terminal. This will remove whatever is in the `Lib` and `Bin` folders.

Thank you for using Eigenious!