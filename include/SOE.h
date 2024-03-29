#ifndef SOE_H
#define SOE_H

#include <vector>
#include "Global_data.h"
#include "Element.h"
using namespace std;

struct SOE
{
    vector<vector<double>> global_H;
    vector<double> global_P;
    vector<vector<double>> global_C;
    vector<double> global_PC;

    SOE(int NodesNumber);
    vector<vector<double>> aggregate_H_C(Element element, Global_data global_data);
    vector<double> aggregate_P(Element element, Global_data global_data);
    vector<double> Results(Global_data global_data, vector<double> InitialTemp);
    vector<double> solveSystem(vector<vector<double>> &global_H, vector<double> &global_P, Global_data global_data);
    void gaussianElimination(vector<vector<double>> &global_H);
};

#endif