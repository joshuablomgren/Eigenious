#include <iostream>
#include "molecule.h"
#include "CNDO.h"

int main () {
    Molecule O2("../Molecules/PS4/O2.txt");
    // Triplet multiplicity
    O2.pAlpha = 7;
    O2.qBeta = 5;

    Molecule O("../Molecules/PS4/O.txt");
    // Triplet multiplicity
    O.pAlpha = 4;
    O.qBeta = 2;

    Matrix<double> S_O2 = calcOverlapMatrix(O2);
    CNDO cndo_O2(O2, S_O2);
    Matrix<double> S_O = calcOverlapMatrix(O);
    CNDO cndo_O(O, S_O);

    cout << "O2: " << endl;
    cout << "pAlpha: " << O2.pAlpha << endl;
    cout << "qBeta: " << O2.qBeta << endl;
    cndo_O2.scfCycle();
    cout << endl;

    cout << "O: " << endl;
    cout << "pAlpha: " << O.pAlpha << endl;
    cout << "qBeta: " << O.qBeta << endl;
    cndo_O.scfCycle();
    cout << endl;

    // Calculate bond energy 2 O -> O2
    cout << "----------------------------------------" << endl;
    double bondEnergy = cndo_O2.totalEnergy - 2 * cndo_O.totalEnergy;
    cout << "Bond Energy 2 O -> O2: " << bondEnergy << " eV" << endl;

    return 0;
}