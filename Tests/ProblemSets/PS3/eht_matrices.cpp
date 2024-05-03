#include "eht_matrices.h"
#include <map>

/**
 * @brief Make the Hamiltonian matrix from the overlap matrix.
 * 
 * @param overlap_mat The overlap matrix
 * 
 * @return The Hamiltonian matrix
 */
Matrix<double> calcHamiltonianMatrix(Molecule molecule, Matrix<double> overlap_mat) {
    // Initialize the Hamiltonian matrix
    Matrix<double> hamiltonian_mat = zeros<double>(molecule.nBasisFunctions, molecule.nBasisFunctions);

    // Define a map to store the atom type to Hamiltonian value mappings
    map<string, double> hamiltonianValues = {
        {"H-1s", -13.6},
        {"C-2s", -21.4},
        {"C-2px", -11.4},
        {"C-2py", -11.4},
        {"C-2pz", -11.4}
    };

    // Loop over all the basis functions of the molecule
    for (int i = 0; i < molecule.nBasisFunctions; i++) {
        for (int j = 0; j < molecule.nBasisFunctions; j++) {
            string atomType1 = molecule.basisFunctionsList[i].AO_type;
            string atomType2 = molecule.basisFunctionsList[j].AO_type;

            // Get the Hamiltonian values given the atom types (i.e H-1s, C-2s, etc.)
            double h_i = hamiltonianValues[atomType1];
            double h_j = hamiltonianValues[atomType2];

            // Assign diagonal elements
            if (i == j) {
                hamiltonian_mat(i, j) = h_i;
            // Calculate off-diagonal elements
            } else {
                hamiltonian_mat(i, j) = constant_K / 2 * (h_i + h_j) * overlap_mat(i, j);
            }
        }
    }
    return hamiltonian_mat;
}

/**
 * @brief Calculate the X transformation matrix by diagonalizing the overlap matrix.
 * 
 * @param overlap_mat The overlap matrix
 * 
 * @return The X transformation matrix
 */
Matrix<double> calcXMatrix(Matrix<double> overlap_mat) {
    // Initialize the X matrix
    Matrix<double> X_mat = zeros<double>(overlap_mat.nRows(), overlap_mat.nCols());

    // Diagonalize the overlap matrix
    Vector<double> s_eigval;
    Matrix<double> U_eigvec;
    eigen_sym(s_eigval, U_eigvec, overlap_mat);

    // Get inverse of square root of eigenvalues
    Vector<double> s_eigval_inv_sqrt = 1.0 / sqrt(s_eigval); 

    // Diagonalize s_eigval_inv_sqrt
    Matrix<double> s_diag_mat = diagonal_matrix(s_eigval_inv_sqrt);

    // Calculate X matrix
    X_mat = U_eigvec * s_diag_mat * U_eigvec.transpose();

    return X_mat;
}

/**
 * @brief Calculate Hamiltonian Prime matrix from the Hamiltonian matrix and the X matrix.
 * 
 * @param X_mat The X transformation matrix
 * @param hamiltonian_mat The Hamiltonian matrix
 * 
 * @return The Hamiltonian Prime matrix
 */
Matrix<double> calcHamiltonianPrimeMatrix(Matrix<double> X_mat, Matrix<double> hamiltonian_mat) {
    // Initialize the Hamiltonian Prime matrix
    Matrix<double> h_prime_mat = zeros<double>(hamiltonian_mat.nRows(), hamiltonian_mat.nCols());

    // Calculate Hamiltonian Prime matrix
    h_prime_mat = X_mat.transpose() * hamiltonian_mat * X_mat;

    return h_prime_mat;
}

/**
 * @brief Calculate the energy of the molecule given the Hamiltonian Prime matrix.
 * 
 * Calculate the energy of the molecule by getting the eigenvalues and eigenvectors of the Hamiltonian Prime matrix.
 * 
 * @param hamiltonian_prime_mat The Hamiltonian Prime matrix
 * 
 * @return The energy of the molecule
 */
double calcEnergy(Matrix<double> X_mat, Matrix<double> hamiltonian_prime_mat, int nElectrons) {
    // Diagonalize the Hamiltonian Prime matrix
    Vector<double> energy_eigval;
    Matrix<double> V_eigvec;
    eigen_sym(energy_eigval, V_eigvec, hamiltonian_prime_mat);

    // Form the MO coffecient matrix
    Matrix<double> C_mat = X_mat * V_eigvec;

    // Calculate the energy by summing over the occupied orbitals
    double energy = 0;
    for (int i = 0; i < nElectrons; i++) {
        energy += 2 * energy_eigval(i);
    }

    return energy;
}

/**
 * @brief Calculate the energy of the molecule given the molecule object and the overlap matrix.
 * 
 * This function wraps the symmetric orthogonalization method to calculate the energy of the molecule.
 * 
 * @param molecule The molecule object
 * @param overlap_mat The overlap matrix
 * 
 * @return The energy of the molecule
 */
double calcEnergy(Molecule molecule, Matrix<double> overlap_mat) {
    Matrix<double> hamiltonian_mat = calcHamiltonianMatrix(molecule, overlap_mat);
    Matrix<double> X_mat = calcXMatrix(overlap_mat);
    Matrix<double> hamiltonian_prime_mat = calcHamiltonianPrimeMatrix(X_mat, hamiltonian_mat);
    double energy = calcEnergy(X_mat, hamiltonian_prime_mat, molecule.nElectrons);

    return energy;
}