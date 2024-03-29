#ifndef SURFACE_H
#define SURFACE_H

#include <vector>
using namespace std;

struct Surface
{
    vector<vector<double>> N;
    Surface(int n)
    {
        N = vector<vector<double>>(n, vector<double>(4));
    }
};

#endif