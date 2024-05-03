#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include "molecule.h"

using namespace std;

const double hartree2eV = 27.2114079527;

class CNDO {
    public:
        CNDO(Molecule molecule, Matrix<double> overlapMatrix);
        Molecule molecule;
        Matrix<double> overlapMatrix;
        Matrix<double> gammaMatrix;
        Matrix<double> hCoreMat;
        Matrix<double> alphaFockMat;
        Matrix<double> betaFockMat;
        Matrix<double> alphaDensityMat;
        Matrix<double> betaDensityMat;
        Matrix<double> oldAlphaDensityMat;
        Matrix<double> oldBetaDensityMat;
        Vector<double> totalDensity;
        Matrix<double> alphaCoeffMat;
        Matrix<double> betaCoeffMat;
        Vector<double> alphaEnergy;
        Vector<double> betaEnergy;
        double nuclearRepulsionEnergy;
        double totalEnergy;

        // Map of semi-empirical parameters 
        map<string, map<string, double>> diagCNDOPara;
        map<string, double> offDiagCNDOPara;

        // Map the AO index to the atom it belongs to
        map<int, int> aoIndexToAtom;

        void scfCycle();

    private:
        int calcNumAOs(string chemSym);

        // Functions used to calculate gamma (2e- integral)
        double calc2eIntegral(AO AO1, AO AO2);
        double pg2eIntegral(Vector<double> center_a, Vector<double> center_b, double sigmaA, double sigmaB);

        Matrix<double> calcGammaMatrix();
        Vector<double> calcTotalDensity(); 
        Matrix<double> calcHCoreMat();
        Matrix<double> calcFockMat(Matrix<double> densityMat);
        Matrix<double> calcDensityMat(Matrix<double> coeffMatA, string type);
        
        double calcNuclearRepulsionEnergy();
        double calcTotalEnergy();
}; 