#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include "molecule.h"

using namespace std;

// Constants
const double constant_K = 1.75;

// Functions used for symmetric orthogonalization to calculate the energy of the molecule
Matrix<double> calcHamiltonianMatrix(Molecule molecule, Matrix<double> overlap_mat);
Matrix<double> calcXMatrix(Matrix<double> overlap_mat);
Matrix<double> calcHamiltonianPrimeMatrix(Matrix<double> X_mat, Matrix<double> hamiltonian_mat);
double calcEnergy(Matrix<double> X_mat, Matrix<double> hamiltonian_prime_mat, int nElectrons);
double calcEnergy(Molecule molecule, Matrix<double> overlap_mat);
