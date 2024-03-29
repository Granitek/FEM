#ifndef ELEMENT_H
#define ELEMENT_H

#include <vector>
#include "MatrixFunctions.h"
#include "Node.h"
using namespace std;

struct Element
{
    int ID[4];
    vector<vector<double>> H;
    vector<vector<double>> HBC;
    vector<vector<double>> C;
    double P[4];

    vector<vector<double>> calculateMatrixH(int n, UniversalElement eu, Node nodes[4], int Conductivity, int SpecificHeat, int Density);
    vector<vector<double>> calculateMatrixHBC(int n, UniversalElement eu, Node nodes[4], int Alfa, int Tot);
};

#endif