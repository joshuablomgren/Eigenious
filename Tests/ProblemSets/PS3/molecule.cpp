#include "molecule.h"
#include "factorial.h"

/**
 * @brief BasisFunction struct constructor
 * 
 */
BasisFunction::BasisFunction(string AO_type, Vector<double> center, Vector<int> lmn, 
                             Vector<double> exponents, Vector<double> contraction_coeffs) {
    this->AO_type = AO_type;
    this->center = center;
    this->lmn = lmn;
    this->exponents = exponents;
    this->contraction_coeffs = contraction_coeffs;

    // Calculate the normalization constants
    calcNormConstants();
}

/**
 * @brief Calculate the three normalization constants for a given basis function.
 * 
 * This function calculates the three normalization constants for a given basis function
 * by evaluating the overlap integral of the basis function with itself. The normalization constants
 * are stored in the basis function struct.
 * 
 */
void BasisFunction::calcNormConstants() {
    // Calculate the overlap integrals for each exponent with itself
    double selfOverlap_1 = overlapIntegral3D(center, center, exponents(0), exponents(0), lmn, lmn);
    double selfOverlap_2 = overlapIntegral3D(center, center, exponents(1), exponents(1), lmn, lmn);
    double selfOverlap_3 = overlapIntegral3D(center, center, exponents(2), exponents(2), lmn, lmn);  

    // Calculate the normalization constants
    norm_constants = {1.0 / sqrt(selfOverlap_1), 1.0 / sqrt(selfOverlap_2), 1.0 / sqrt(selfOverlap_3)};
}

/**
 * @brief Molecule class constructor
 * 
 * This constructor reads in the atomic information from a file, and sets up the basis
 * functions and number of electrons of the molecule.
 * 
 * @param filename The name of the file to read from
 * 
 */
Molecule::Molecule(const char* filename) {
    // Open the file
    ifstream infile(filename);
    assert(infile.good());

    // Read in the number of atoms
    infile >> nAtoms;

    atomicNumbers.resize(nAtoms);
    coordinates.resize(nAtoms, 3);

    std::string line; 
    getline(infile, line); // Skip the rest of the line

    for (int i = 0; i < nAtoms; i++) {
        // Get line and split it
        getline(infile, line);
        istringstream iss(line);
        iss >> atomicNumbers(i);
        // Loop over the three coordinates x, y, z
        for (int j = 0; j < 3; j++) {
            iss >> coordinates(i, j);
        }
    }

    // Close the file
    infile.close();

    // Count carbon and hydrogen atoms
    carbon_count = 0;
    hydrogen_count = 0;
    for (int i = 0; i < nAtoms; i++) {
        if (atomicNumbers(i) == 1) {
            hydrogen_count++;
        } else if (atomicNumbers(i) == 6) {
            carbon_count++;
        }
    }

    // Calculate the number of electrons and basis functions
    nElectrons = countElectrons();
    nBasisFunctions = countBasisFunctions();

    // Information on H and C basis functions
    H_exponents = {3.42525091, 0.62391373, 0.16885540};
    H_coefficients = {0.15432897, 0.53532814, 0.44463454};
    C_exponents = {2.94124940, 0.68348310, 0.22228990};
    C_2s_coefficients = {-0.09996723, 0.39951283, 0.70011547};
    C_2p_coefficients = {0.15591627, 0.60768372, 0.39195739};

    // Build the list of basic functions
    basisFunctionsList = buildBasisFunctionsList();
}

/**
 * @brief Evaluate the number of basis functions in the molecule
 * 
 * @return int The number of electrons
 */
int Molecule::countBasisFunctions() {
    return 4 * carbon_count + hydrogen_count;
}

/**
 * @brief Evaluate the number of electrons in the molecule
 * 
 * @return int The number of electrons
 */
int Molecule::countElectrons() {
    double electron_pairs = 2 * carbon_count + (hydrogen_count / 2);

    // Check if electron_pairs is not an integer
    if (electron_pairs != int(electron_pairs)) {
        cout << "Error: Number of electron pairs is not an integer." << endl;
        exit(1);
    } else {
        return int(electron_pairs);
    }
}


/**
 * @brief Create a list of all the basis functions in the molecule, which are contracted Gaussians.
 * 
 * This function creates a list of all the basis functions in the molecule, which are contracted Gaussians.
 * Each basis function is stored as a struct, which contains the type of atomic orbital, the center, the
 * angular momentum, the exponents, the contraction coefficients, and the normalization constants.
 * 
 * @return vector<BasisFunction> A vector of basic functions
 */
vector<BasisFunction> Molecule::buildBasisFunctionsList() {
    // Initialize the list of basis functions
    vector<BasisFunction> basisFunctionsList;

    // Loop over all atoms
    for (int i = 0; i < nAtoms; i++) {

        // If the atom is hydrogen, add a 1s basis function
        if (atomicNumbers(i) == 1) {
            struct BasisFunction AO_1s = {"H-1s", coordinates.getRow(i), {0, 0, 0}, H_exponents, H_coefficients};
            basisFunctionsList.push_back(AO_1s);

        // If the atom is carbon, add 2s and 2p basis functions (4 in total)
        } else if (atomicNumbers(i) == 6) {
            struct BasisFunction AO_2s = {"C-2s", coordinates.getRow(i), {0, 0, 0}, C_exponents, C_2s_coefficients};
            struct BasisFunction AO_2px = {"C-2px", coordinates.getRow(i), {1, 0, 0}, C_exponents, C_2p_coefficients};
            struct BasisFunction AO_2py = {"C-2py", coordinates.getRow(i), {0, 1, 0}, C_exponents, C_2p_coefficients};
            struct BasisFunction AO_2pz = {"C-2pz", coordinates.getRow(i), {0, 0, 1}, C_exponents, C_2p_coefficients};
            basisFunctionsList.push_back(AO_2s);
            basisFunctionsList.push_back(AO_2px);
            basisFunctionsList.push_back(AO_2py);
            basisFunctionsList.push_back(AO_2pz);
        }
    }
    return basisFunctionsList;
}

/**
 * @brief Print the molecule information
 * 
 */
void Molecule::printMoleculeInfo() {
    cout << "Number of atoms: " << nAtoms << endl;
    cout << "Atomic numbers: " << atomicNumbers << endl;
    cout << "Coordinates:\n" << coordinates << endl;
    cout << "Carbon count: " << carbon_count << endl;
    cout << "Hydrogen count: " << hydrogen_count << endl;
    cout << "Number of electrons: " << nElectrons << endl;
    cout << "Number of basis functions: " << nBasisFunctions << endl;
}


/**
 * @brief Calculate the analytical integration for the overlap of two primitive Gaussian shells for one dimension.
 * 
 * This function calculates the analytical integration for the overlap of two primitive Gaussian shells
 * for one dimension given the centers, alpha exponents, and angular momentum of the two shells.
 * 
 * @param alpha The exponent of the first primitive Gaussian shell
 * @param beta The exponent of the second primitive Gaussian shell
 * @param center_a The center of the first primitive Gaussian shell for one dimension
 * @param center_b The center of the second primitive Gaussian shell for one dimension
 * @param lA The angular momentum of the first primitive Gaussian shell
 * @param lB The angular momentum of the second primitive Gaussian shell
 * 
 * @return double The analytical overlap integral
 */
double overlapIntegral1D(double alpha, double beta, double center_a, double center_b, int lA, int lB) {
    // Calculate the exponential prefactor and the associated square root 
    double prefactor = exp(-alpha * beta * pow(center_a - center_b, 2) / (alpha + beta));
    prefactor *= sqrt(M_PI / (alpha + beta));

    // Calculate center_product
    double center_product = (alpha * center_a + beta * center_b) / (alpha + beta);

    // Double summation over the angular momentum combinations
    double sum = 0.0;
    for (int i = 0; i <= lA; i++) {
        for (int j = 0; j <= lB; j++) {
            // Only (i + j) even terms contribute
            if ((i + j) % 2 == 0) {
                sum += binomialCoef(lA, i) * binomialCoef(lB, j) * (doubleFactorial(i + j - 1) 
                        * pow(center_product - center_a, lA - i) * pow(center_product - center_b, lB - j)) 
                        / pow(2 * (alpha + beta), double(i + j) / 2);
            }
        }
    }
    double integral = prefactor * sum;
    return integral;
}

/**
 * @brief Calculate the 3D overlap integral for two primitive Gaussian shells.
 * 
 * This function evaluates the overlap integral of two primitive Gaussian shells over all dimensions and 
 * over all angular momentum combinations.
 * 
 * @param centers_a The center of the first primitive Gaussian shell for all dimensions
 * @param centers_b The center of the second primitive Gaussian shell for all dimensions
 * @param alpha The exponent of the first primitive Gaussian shell
 * @param beta The exponent of the second primitive Gaussian shell
 * @param lmn_a The angular momentum of the first primitive Gaussian shell
 * @param lmn_b The angular momentum of the second primitive Gaussian shell
 * 
 * @return double The overlap integral 
 */
double overlapIntegral3D(Vector<double> centers_a, Vector<double> centers_b, double alpha, double beta, Vector<int> lmn_a, Vector<int> lmn_b) {
    // Calculate the overlap integral for each dimension
    double integral = overlapIntegral1D(alpha, beta, centers_a(0), centers_b(0), lmn_a(0), lmn_b(0)) *
                      overlapIntegral1D(alpha, beta, centers_a(1), centers_b(1), lmn_a(1), lmn_b(1)) *
                      overlapIntegral1D(alpha, beta, centers_a(2), centers_b(2), lmn_a(2), lmn_b(2));

    return integral;
}

/**
 * @brief Calculate the contracted overlap integral for two basis functions.
 * 
 * This function calculates the contracted overlap integral for two basis functions after the
 * normalization constants have been calculated.
 * 
 * @param basisFunction1 The first basis function
 * @param basisFunction2 The second basis function
 * 
 * @return double The contracted overlap integral
 */
double calcContractedOverlap(BasisFunction basisFunction1, BasisFunction basisFunction2) {
    double contracted_overlap = 0.0;

    // Loop over all exponents
    for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
            double unnorm_overlap = overlapIntegral3D(basisFunction1.center, basisFunction2.center, 
                                                      basisFunction1.exponents(k), basisFunction2.exponents(l), 
                                                      basisFunction1.lmn, basisFunction2.lmn);
                                                      
            contracted_overlap += basisFunction1.contraction_coeffs(k) * basisFunction2.contraction_coeffs(l) *
                                  basisFunction1.norm_constants(k) * basisFunction2.norm_constants(l) * unnorm_overlap;
        }
    }
    return contracted_overlap;
}

/**s
 * @brief Calculate the contract overlap matrix for the entire molecule.
 * 
 * This function calculates the contracted overlap matrix for the entire molecule given a molecule by
 * looping over all basis function combinations. 
 * 
 * @param molecule The molecule to calculate the overlap matrix for
 * 
 * @return Matrix<double> The overlap matrix
 */
Matrix<double> calcOverlapMatrix(Molecule molecule) {
    // Initialize the overlap matrix
    Matrix<double> overlap_matrix = zeros<double>(molecule.nBasisFunctions, molecule.nBasisFunctions);

    // Loop over all basis function combinations
    for (int i = 0; i < molecule.nBasisFunctions; i++) {
        for (int j = 0; j < molecule.nBasisFunctions; j++) {
            overlap_matrix(i, j) = calcContractedOverlap(molecule.basisFunctionsList[i], molecule.basisFunctionsList[j]);
        }
    }

    return overlap_matrix;
}



