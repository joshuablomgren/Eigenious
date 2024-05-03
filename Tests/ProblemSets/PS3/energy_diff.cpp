#include <iostream>
#include "molecule.h"
#include "eht_matrices.h"

// Evaluating chemical energy differences

int main () {

    // Finding the bond energy of H2
    cout << "Finding Bond Energy of H2: " << endl << endl;
    Molecule H2("../Molecules/PS3/H2.txt");

    // Matrices and Energy of H2
    Matrix<double> overlap_mat_H2 = calcOverlapMatrix(H2);
    Matrix<double> hamiltonian_mat_H2 = calcHamiltonianMatrix(H2, overlap_mat_H2);
    double energy_H2 = calcEnergy(H2, overlap_mat_H2);
    cout << "Energy of H2: " << energy_H2 << "eV" << endl;

    // Subtract twice the energy of the H atom
    double energy_H = -13.6;
    double H2_bond_energy = energy_H2 - 2 * energy_H;

    // Print energy difference
    cout << "Bond energy of H2: " << H2_bond_energy << "eV" << endl;


    cout << "----------------------------------------" << endl;
    cout << "Finding Energy Difference: C2H2 + H2 <-> C2H4" << endl << endl;


    // Finding the energy difference for the chemical reaction C2H2 + H2 -> C2H4
    Molecule C2H2("../Molecules/PS3/C2H2.txt");
    Molecule C2H4("../Molecules/PS3/C2H4.txt");
    Molecule H2_2("../Molecules/PS3/H2.txt");

    // Energies of C2H2, C2H4, and H2
    Matrix<double> overlap_mat_C2H2 = calcOverlapMatrix(C2H2);
    double energy_C2H2 = calcEnergy(C2H2, overlap_mat_C2H2);

    Matrix<double> overlap_mat_C2H4 = calcOverlapMatrix(C2H4);
    double energy_C2H4 = calcEnergy(C2H4, overlap_mat_C2H4);

    Matrix<double> overlap_mat_H2_2 = calcOverlapMatrix(H2_2);
    double energy_H2_2 = calcEnergy(H2_2, overlap_mat_H2_2);

    // Print energies
    cout << "Energy of C2H4: " << energy_C2H4 << "eV" << endl;
    cout << "Energy of C2H2: " << energy_C2H2 << "eV" << endl;
    cout << "Energy of H2: " << energy_H2_2 << "eV" << endl << endl;

    // Subtract the energies of the reactants from the energy of the product
    double C2H4_energy_diff = energy_C2H4 - energy_C2H2 - energy_H2_2;

    // Convert to kJ/mol
    C2H4_energy_diff *= 96.485;

    // Print energy difference
    cout << "Energy difference for the chemical reaction C2H2 + H2 <-> C2H4: " << C2H4_energy_diff << "eV" << endl;

    return 0;
}