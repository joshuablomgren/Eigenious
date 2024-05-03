#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include "../../../Include/matrix.h"
#include "../../../Include/vector.h"

using namespace std;

struct AO {
    string AO_type;
    string chemSym;
    int valence;
    Vector<double> center;
    Vector<int> lmn;
    Vector<double> exponents;
    Vector<double> contraction_coeffs;
    Vector<double> norm_constants;

    AO(string AO_type, string chemSym, int valence, Vector<double> center, \
       Vector<int> lmn, Vector<double> exponents, Vector<double> contraction_coeffs);
    void calcNormConstants();
};

class Molecule {
    public:
        Molecule(const char* filename);
        int nAtoms;
        Vector<int> atomicNumbers;
        vector<string> atomicSymbols;
        Vector<int> atomValences;
        Matrix<double> coordinates; 
        int charge;
        
        int hydrogenCount;
        int heavyAtomCount;
        int nElectrons; 
        int pAlpha;
        int qBeta;
        int nBasisFunctions;
        vector<AO> basisFunctionsList;
        map<int, int> aoIndexToAtom;

        map<string, Vector<double>> exponentsMap;
        map<string, Vector<double>> sContractionMap;
        map<string, Vector<double>> pContractionMap;

        void printMoleculeInfo();

    private:
        int countBasisFunctions();
        int countElectrons();
        vector<AO> buildBasisFunctionsList();
};

// Functions used to calculate overlap integrals
double overlapIntegral1D(double alpha, double beta, double center_a, double center_b, int lA, int lB);
double overlapIntegral3D(Vector<double> centers_a, Vector<double> centers_b, double alpha, double beta, Vector<int> lmn_a, Vector<int> lmn_b);
double calcContractedOverlap(AO AO1, AO AO2);
Matrix<double> calcOverlapMatrix(Molecule molecule);
