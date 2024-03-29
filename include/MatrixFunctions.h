#ifndef MATRIXFUNCTIONS_H
#define MATRIXFUNCTIONS_H

#include <vector>
#include <iostream>
#include "UniversalElement.h"
using namespace std;

void printMatrix(const vector<vector<double>> &matrix);
vector<vector<double>> transposeMatrix(const vector<vector<double>> &matrix);
vector<vector<vector<double>>> transposeSurfaces(const UniversalElement &eu);

#endif