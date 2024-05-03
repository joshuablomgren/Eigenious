// main.cpp
#include <iostream>
#include "molecule.h"
#include "eht_matrices.h"
#include <chrono>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>\n";
        return 1;
    }
    // Time the program
    auto start_time = chrono::high_resolution_clock::now();

    const char* filename = argv[1];

    Molecule molecule(filename);

    molecule.printMoleculeInfo();

    cout << "----------------------------------------" << endl << endl;
    
    // Test contracted overlap integral
    Matrix<double> overlap_mat = calcOverlapMatrix(molecule);
    cout << "Overlap matrix:\n" << overlap_mat << endl;

    // Test Hamiltonian matrix
    Matrix<double> hamiltonian_mat = calcHamiltonianMatrix(molecule, overlap_mat);
    cout << "Hamiltonian matrix:\n" << hamiltonian_mat << endl;

    // Test X matrix
    Matrix<double> X_mat = calcXMatrix(overlap_mat);
    cout << "X matrix:\n" << X_mat << endl;

    // Test Hamiltonian Prime matrix
    Matrix<double> hamiltonian_prime_mat = calcHamiltonianPrimeMatrix(X_mat, hamiltonian_mat);
    cout << "Hamiltonian Prime matrix:\n" << hamiltonian_prime_mat << endl;

    // Test calcEnergy
    double energy = calcEnergy(molecule, overlap_mat);
    cout << "Energy: " << energy << endl;

    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << "Time taken: " << duration.count() << " microseconds" << endl;

    return 0;
}
