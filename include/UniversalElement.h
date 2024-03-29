#ifndef UNIVERSALELEMENT_H
#define UNIVERSALELEMENT_H

#include <vector>
#include <functional>
#include <iostream>
#include "Surface.h"
#include "GausTable.h"

using namespace std;

struct UniversalElement
{
    vector<vector<double>> ksi;
    vector<vector<double>> eta;
    vector<Surface> Surface3D;
    vector<vector<double>> N_nodes;
    vector<function<double(double)>> ksi_function = vector<function<double(double)>>();
    vector<function<double(double)>> eta_function = vector<function<double(double)>>();
    vector<function<double(double, double)>> N_function = vector<function<double(double, double)>>();

    UniversalElement(int n);
    void shape_functions();
    void N_data(int n);
    void N_data_nodes(int n);
    void print_Surface3D(int n);
    void ksi_eta_data(int n);
    void print(int n);
};

#endif