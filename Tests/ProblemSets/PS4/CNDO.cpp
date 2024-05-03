#include "CNDO.h"

/**
 * @brief CNDO/2 Hartree Fock method constructor
 * 
 */
CNDO::CNDO(Molecule molecule, Matrix<double> overlapMatrix)
: molecule(molecule), overlapMatrix(overlapMatrix) {

    // Map the AO index to the atom it belongs to
    int index = 0;
    for (int A = 0; A < molecule.nAtoms; A++) {
        int numAOs = calcNumAOs(molecule.atomicSymbols[A]);
        for (int i = 0; i < numAOs; i++) {
            aoIndexToAtom[index] = A;
            index++;
        }
    }

    // Initialize semi-empirical parameters
    diagCNDOPara["H"]["1s"] = 7.176;

    diagCNDOPara["C"]["2s"] = 14.051;
    diagCNDOPara["C"]["2px"] = 5.572;
    diagCNDOPara["C"]["2py"] = 5.572;
    diagCNDOPara["C"]["2pz"] = 5.572;

    diagCNDOPara["N"]["2s"] = 19.316;
    diagCNDOPara["N"]["2px"] = 7.275;
    diagCNDOPara["N"]["2py"] = 7.275;
    diagCNDOPara["N"]["2pz"] = 7.275;

    diagCNDOPara["O"]["2s"] = 25.390;
    diagCNDOPara["O"]["2px"] = 9.111;
    diagCNDOPara["O"]["2py"] = 9.111;
    diagCNDOPara["O"]["2pz"] = 9.111;

    diagCNDOPara["F"]["2s"] = 32.272;
    diagCNDOPara["F"]["2px"] = 11.080;
    diagCNDOPara["F"]["2py"] = 11.080;
    diagCNDOPara["F"]["2pz"] = 11.080;

    offDiagCNDOPara["H"] = 9.;
    offDiagCNDOPara["C"] = 21.;
    offDiagCNDOPara["N"] = 25.;
    offDiagCNDOPara["O"] = 31.;
    offDiagCNDOPara["F"] = 39.;

    // Initialize matrices
    alphaCoeffMat = zeros<double>(molecule.nBasisFunctions, molecule.nBasisFunctions);
    betaCoeffMat = zeros<double>(molecule.nBasisFunctions, molecule.nBasisFunctions);

    alphaDensityMat = zeros<double>(molecule.nBasisFunctions, molecule.nBasisFunctions);
    betaDensityMat = zeros<double>(molecule.nBasisFunctions, molecule.nBasisFunctions);
    totalDensity = zeros<double>(molecule.nBasisFunctions);

    alphaEnergy.reset_size(molecule.nBasisFunctions);
    betaEnergy.reset_size(molecule.nBasisFunctions);

    gammaMatrix = calcGammaMatrix();
    alphaFockMat = calcFockMat(alphaDensityMat);
    betaFockMat = calcFockMat(betaDensityMat);
    hCoreMat = calcHCoreMat();

    nuclearRepulsionEnergy = calcNuclearRepulsionEnergy();
}

/**
 * @brief Calculate the number of AOs for an atom.
 * 
 * This function calculates the number of AOs for an atom given the atomic
 * symbol. Hydrogen atoms have 1 AO, while all other atoms have 4 AOs.
 * 
 * @param chemSym The atomic symbol of the atom
 * 
 * @return int The number of AOs
 */
int CNDO::calcNumAOs(string chemSym) {
    if (chemSym == "H") {
        return 1;
    }
    else {
        return 4;
    }
}

/**
 * @brief Calculate the 2-electron repulsion integrals for two primitive Gaussian shells.
 * 
 * This function calculates the 2-electron repulsion integrals evaluated over the
 * square of the valence s orbital centered on atoms A and B.
 * 
 * @param AO1 The first AO
 * @param AO2 The second AO
 * 
 * @return double The 2-electron repulsion integral
 */
double CNDO::calc2eIntegral(AO AO1, AO AO2) {
    if (!(AO1.lmn.accumulate() == 0 && AO2.lmn.accumulate() == 0)) {
        cout << "Error: 2e integrals only implemented for s orbitals" << endl;
    }

    // d prime = contraction coefficients * normalization constant
    Vector<double> dprime_a = AO1.contraction_coeffs % AO1.norm_constants;
    Vector<double> dprime_b = AO2.contraction_coeffs % AO2.norm_constants;

    int len = AO1.exponents.size();

    double gamma = 0.0;
    for (int k1 = 0; k1 < len; k1++) {
        for (int k2 = 0; k2 < len; k2++) {
            double sigmaA = 1.0/(AO1.exponents(k1) + AO1.exponents(k2)); // eq 3.10
            for (int l1 = 0; l1 < len; l1++) {
                for (int l2 = 0; l2 < len; l2++) {
                    double sigmaB = 1.0/(AO2.exponents(l1) + AO2.exponents(l2)); // eq 3.10
                    
                    double I2e = pg2eIntegral(AO1.center, AO2.center, sigmaA, sigmaB);  // eq 3.14
                    gamma += dprime_a(k1) * dprime_a(k2) * dprime_b(l1) * dprime_b(l2) * I2e;
                }
            }
        }
    }
    return gamma;
}

/**
 * @brief Calculates the 2-electron integral over primitive Gaussians.
 * 
 * This function calculates the 2-electron integral over primitive Gaussians
 * given the centers, exponents, and contraction coefficients of the two Gaussians.
 * 
 * @param center_a 
 * @param center_b 
 * @param sigmaA 
 * @param sigmaB 
 * 
 * @return double 2-electron integral
 */
double CNDO::pg2eIntegral(Vector<double> center_a, Vector<double> center_b, double sigmaA, double sigmaB) {
    double U = pow(M_PI * sigmaA, 1.5) * pow(M_PI * sigmaB, 1.5);
    double V2 = 1.0 / (sigmaA + sigmaB);

    double distance = norm(center_a - center_b);

    if (distance == 0.0) {
        return U * sqrt(2 * V2) * sqrt(2 / M_PI) * hartree2eV;  // eq 3.15
    } 

    double sqrtT = sqrt(V2) * distance;
    double result = U / distance * erf(sqrtT); 
    return result * hartree2eV;
}

/**
 * @brief Calculate the gamma matrix of 2-electron repulsion integrals.
 * 
 * This function calculates the gamma matrix of 2-electron repulsion integrals
 * for the entire molecule using the calc2eIntegral function.
 * 
 * @param molecule 
 * 
 * @return mat The gamma matrix
 */
Matrix<double> CNDO::calcGammaMatrix() {
    // Make a list of only basis functions with s orbitals
    vector<AO> sBasisFunctionsList;
    for (int i = 0; i < molecule.nBasisFunctions; i++) {
        if (molecule.basisFunctionsList[i].AO_type == "1s" || molecule.basisFunctionsList[i].AO_type == "2s") {
            sBasisFunctionsList.push_back(molecule.basisFunctionsList[i]);
        }
    }

    Matrix<double> gamma_matrix = zeros<double>(molecule.nAtoms, molecule.nAtoms);

    // Loop over all s orbital basis function combinations
    for (int i = 0; i < sBasisFunctionsList.size(); i++) {
        for (int j = 0; j < sBasisFunctionsList.size(); j++) {
            gamma_matrix(i, j) = calc2eIntegral(sBasisFunctionsList[i], sBasisFunctionsList[j]);
        }
    }
    return gamma_matrix;
}

/**
 * @brief Calculate the nuclear repulsion energy of the molecule.
 * 
 * This function calculates the nuclear repulsion energy of the molecule
 * given the charges of the atoms
 * 
 * @return double The nuclear repulsion energy
 */
double CNDO::calcNuclearRepulsionEnergy() {
    double nuclearRepulsionEnergy = 0.0;
    for (int A = 0; A < molecule.nAtoms; A++) {
        for (int B = 0; B < A; B++) {
            double distance = norm(molecule.coordinates.getRow(A) - molecule.coordinates.getRow(B));
            nuclearRepulsionEnergy += molecule.atomValences(A) * molecule.atomValences(B) / distance;
        }
    }
    return nuclearRepulsionEnergy * hartree2eV;
}

   
/**
 * @brief Calculate the total density vector.
 * 
 * This function calculates the total density vector using the alpha and beta
 * density matrices.
 * 
 * @return mat The total density vector
 */
Vector<double> CNDO::calcTotalDensity() {
    Vector<double> totalDensity = zeros<double>(molecule.nAtoms);

    for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
        // Get atom associated with mu AO
        int A = aoIndexToAtom[mu];
        totalDensity(A) += alphaDensityMat(mu, mu) + betaDensityMat(mu, mu);
    }

    return totalDensity;
}

/**
 * @brief Calculate the CNDO/2 Fock matrix.
 * 
 * This function calculates the CNDO/2 Fock matrix given either the alpha or
 * beta density matrix.
 * 
 * @param densityMat The density matrix (alpha or beta)
 * 
 * @return mat The CNDO/2 Fock matrix
 */
Matrix<double> CNDO::calcFockMat(Matrix<double> densityMat) {
    Matrix<double> fockMat = zeros<double>(molecule.nBasisFunctions, molecule.nBasisFunctions);

    // Loop over all AOs in molecule
    for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
        for (int nu = 0; nu < molecule.nBasisFunctions; nu++) {
            // Get atoms associated with mu and nu AOs
            int A = aoIndexToAtom[mu];
            int B = aoIndexToAtom[nu];
            string chemSymA = molecule.atomicSymbols[A];
            string chemSymB = molecule.atomicSymbols[B];
            double gammaAA = gammaMatrix(A, A);
            double gammaAB = gammaMatrix(A, B);
            double pAA = totalDensity(A);
            double ZA = molecule.atomValences(A);

            // Calculate the diagonal elements of the matrix
            if (mu == nu) {
                string AO_type = molecule.basisFunctionsList[mu].AO_type;
                fockMat(mu, nu) = -diagCNDOPara[chemSymA][AO_type] + \
                                  ((pAA - ZA) - (densityMat(mu, mu) - 0.5)) * gammaAA;
                
                // Update the diagonal elements of the matrix when A != B
                for (int B = 0; B < molecule.nAtoms; B++) {
                    if (A != B) {
                        double pBB = totalDensity(B);
                        double ZB = molecule.atomValences(B);
                        double gammaAB = gammaMatrix(A, B);
                        fockMat(mu, nu) += (pBB - ZB) * gammaAB;
                    }
                }
            }

            // Calculate the off-diagonal elements of the matrix
            else {
                fockMat(mu, nu) = (-offDiagCNDOPara[chemSymA] - offDiagCNDOPara[chemSymB]) \
                                  / 2.0 * overlapMatrix(mu, nu) - (densityMat(mu, nu) * gammaAB);
            }
        }
    }
    return fockMat;
}

/**
 * @brief Calculate the core Hamiltonian matrix.
 * 
 * This function calculates the core Hamiltonian matrix similar to the Fock
 * matrix, but independent of electron density.
 * 
 * @return mat The core Hamiltonian matrix
 */
Matrix<double> CNDO::calcHCoreMat() {
    Matrix<double> hCoreMat = zeros<double>(molecule.nBasisFunctions, molecule.nBasisFunctions);
    
    // Loop over all AOs in molecule
    for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
        for (int nu = 0; nu < molecule.nBasisFunctions; nu++) {
            int A = aoIndexToAtom[mu];
            int B = aoIndexToAtom[nu];
            string chemSymA = molecule.atomicSymbols[A];
            string chemSymB = molecule.atomicSymbols[B];
            double gammaAA = gammaMatrix(A, A);
            double gammaAB = gammaMatrix(A, B);
            double ZA = molecule.atomValences(A);

            // Calculate the diagonal elements of the matrix
            if (mu == nu) {
                string AO_type = molecule.basisFunctionsList[mu].AO_type;
                hCoreMat(mu, nu) = -diagCNDOPara[chemSymA][AO_type] - (ZA - 0.5) * gammaAA;

                for (int B = 0; B < molecule.nAtoms; B++) {
                    if (A != B) {
                        double ZB = molecule.atomValences(B);
                        double gammaAB = gammaMatrix(A, B);
                        hCoreMat(mu, nu) -= ZB * gammaAB;
                    }
                }
            }
            // Calculate the off-diagonal elements of the matrix
            else {
                hCoreMat(mu, nu) = (-offDiagCNDOPara[chemSymA] - offDiagCNDOPara[chemSymB]) \
                                   / 2.0 * overlapMatrix(mu, nu);
            }
        }
    }
    return hCoreMat;
}

/**
 * @brief Calculate the density matrix from the coefficient matrix.
 * 
 * This function calculates the density matrix given the coefficient matrix
 * for either the alpha or beta electrons.
 * 
 * @param coeffMat The coefficient matrix to use (alpha or beta)
 * @param type The type of matrix (alpha or beta)
 * 
 * @return mat The density matrix
 */
Matrix<double> CNDO::calcDensityMat(Matrix<double> coeffMatA, string type) {
    Matrix<double> densityMat = zeros<double>(molecule.nBasisFunctions, molecule.nBasisFunctions);

    if (type == "alpha") {
        for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
            for (int nu = 0; nu < molecule.nBasisFunctions; nu++) {
                for (int i = 0; i < molecule.pAlpha; i++) {
                    densityMat(mu, nu) += coeffMatA(mu, i) * coeffMatA(nu, i);
                }
            }
        }
    }

    else if (type == "beta") {
        for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
            for (int nu = 0; nu < molecule.nBasisFunctions; nu++) {
                for (int i = 0; i < molecule.qBeta; i++) {
                    densityMat(mu, nu) += coeffMatA(mu, i) * coeffMatA(nu, i);
                }
            }
        }
    }
    return densityMat;
}

/**
 * @brief Calculate the total energy.
 * 
 * This function calculates the total energy after convergence is reached.
 * 
 * @return double The total energy
 */
double CNDO::calcTotalEnergy() {
    double totalEnergy = 0.0;
    for (int mu = 0; mu < molecule.nBasisFunctions; mu++) {
        for (int nu = 0; nu < molecule.nBasisFunctions; nu++) {
            totalEnergy += (alphaDensityMat(mu, nu) * (hCoreMat(mu, nu) + alphaFockMat(mu, nu)) + \
                            betaDensityMat(mu, nu) * (hCoreMat(mu, nu) + betaFockMat(mu, nu)));
        }
    }
    totalEnergy /= 2.0;
    totalEnergy += nuclearRepulsionEnergy;
    return totalEnergy;
}

/**
 * @brief Self consistent field (SCF) cycle.
 * 
 * This function runs the SCF cycle until convergence is reached. The SCF cycle
 * consists of the following steps:
 * 1) Guess the density matrix = 0
 * 2) Calculate the Fock matrix
 * 3) Diagonalize the Fock matrix and obtain the coefficient matrix
 * 4) Calculate the density matrix
 * 5) Check for convergence
 * 6) If not converged, repeat from step 2
 * 7) If converged, calculate the total energy
 */
void CNDO::scfCycle() {
    // 1) Guess the density matrix
    alphaDensityMat = zeros<double>(molecule.nBasisFunctions, molecule.nBasisFunctions);
    betaDensityMat = zeros<double>(molecule.nBasisFunctions, molecule.nBasisFunctions);

    bool converged = false;
    int scfCycleCount = 0;  

    // 6) If not converged, repeat from step 2
    while (!converged) {
        scfCycleCount++;
        cout << "---------------------------------------------------------" << endl;
        cout << "SCF cycle: " << scfCycleCount << endl << endl;

        // 2) Calculate the Fock matrix
        alphaFockMat = calcFockMat(alphaDensityMat);
        betaFockMat = calcFockMat(betaDensityMat);

        // 3) Diagonalize the Fock matrix and obtain the coefficient matrix
        eigen_sym(alphaEnergy, alphaCoeffMat, alphaFockMat);
        eigen_sym(betaEnergy, betaCoeffMat, betaFockMat);

        // Save
        oldAlphaDensityMat = alphaDensityMat;
        oldBetaDensityMat = betaDensityMat;

        // 4) Calculate the density matrix
        alphaDensityMat = calcDensityMat(alphaCoeffMat, "alpha");
        betaDensityMat = calcDensityMat(betaCoeffMat, "beta");
        totalDensity = calcTotalDensity();

         // Print matrices
        cout << "Alpha Density Mat: " << endl;
        cout << alphaDensityMat << endl;
        cout << "Beta Density Mat: " << endl;
        cout << betaDensityMat << endl;
        cout << "Total Density: " << endl;
        cout << totalDensity << endl;
        cout << "Alpha Fock Mat: " << endl;
        cout << alphaFockMat << endl;
        cout << "Beta Fock Mat: " << endl;
        cout << betaFockMat << endl;
        cout << "Alpha Coeff Mat: " << endl;
        cout << alphaCoeffMat << endl;
        cout << "Beta Coeff Mat: " << endl;
        cout << betaCoeffMat << endl;

        // 5) Check for convergence (tolerance = 1e-6)
        if (abs(alphaDensityMat - oldAlphaDensityMat).max() < 1e-6 && \
            abs(betaDensityMat - oldBetaDensityMat).max() < 1e-6) {
                
            // 7) If converged, calculate the total energy
            converged = true;
            cout << "*********************************************************" << endl;
            cout << "SCF cycle converged after " << scfCycleCount << " iterations!" << endl;
            totalEnergy = calcTotalEnergy();
            cout << "Nuclear Repulsion Energy: " << nuclearRepulsionEnergy << " eV" << endl;
            cout << "Total Energy: " << totalEnergy << " eV" << endl;
        }
    }
}
