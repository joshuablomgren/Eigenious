#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include "../../../Include/matrix.h"
#include "../../../Include/vector.h"

using namespace std;

struct BasisFunction {
    string AO_type;
    Vector<double> center;
    Vector<int> lmn;
    Vector<double> exponents;
    Vector<double> contraction_coeffs;
    Vector<double> norm_constants;

    BasisFunction(string AO_type, Vector<double> center, Vector<int> lmn, \
                  Vector<double> exponents, Vector<double> contraction_coeffs);
    void calcNormConstants();
};

// Class that contains all the information about a molecule
class Molecule {
    public:
        Molecule(const char* filename);
        int nAtoms;
        Vector<int> atomicNumbers;
        Matrix<double> coordinates;
        int carbon_count;
        int hydrogen_count;
        int nElectrons; 
        int nBasisFunctions;
        vector<BasisFunction> basisFunctionsList;

        // Information on H and C basis functions
        Vector<double> H_exponents;
        Vector<double> H_coefficients;
        Vector<double> C_exponents;
        Vector<double> C_2s_coefficients;
        Vector<double> C_2p_coefficients;

        void printMoleculeInfo();

    private:
        int countBasisFunctions();
        int countElectrons();
        vector<BasisFunction> buildBasisFunctionsList();
};

// Functions used to calculate overlap integrals
double overlapIntegral1D(double alpha, double beta, double center_a, double center_b, int lA, int lB);
double overlapIntegral3D(Vector<double> centers_a, Vector<double> centers_b, double alpha, double beta, Vector<int> lmn_a, Vector<int> lmn_b);
double calcContractedOverlap(BasisFunction basisFunction1, BasisFunction basisFunction2);
Matrix<double> calcOverlapMatrix(Molecule molecule);