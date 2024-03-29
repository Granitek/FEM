#ifndef GAUSTABLE_H
#define GAUSTABLE_H

#include <vector>
#include <cmath>
#include <iostream>
using namespace std;

struct GausTable
{
    vector<double> nodes;
    vector<double> weights;

    GausTable(int n);
    double Gauss_solution(int n, int m);
};

#endif