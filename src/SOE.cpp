#include "../include/SOE.h"

SOE::SOE(int NodesNumber)
{
    global_H.resize(NodesNumber, vector<double>(NodesNumber, 0.0));
    global_P.resize(NodesNumber, 0.0);
    global_C.resize(NodesNumber, vector<double>(NodesNumber, 0.0));
    global_PC.resize(NodesNumber, 0.0);
}

vector<vector<double>> SOE::aggregate_H_C(Element element, Global_data global_data)
{
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
        {
            global_H[element.ID[i] - 1][element.ID[j] - 1] += element.H[i][j] + element.HBC[i][j];
            global_C[element.ID[i] - 1][element.ID[j] - 1] += element.C[i][j];
        }
    return global_H;
}

vector<double> SOE::aggregate_P(Element element, Global_data global_data)
{
    for (int i = 0; i < 4; i++)
        global_P[element.ID[i] - 1] += element.P[i];
    return global_P;
}

vector<double> SOE::Results(Global_data global_data, vector<double> InitialTemp)
{
    vector<vector<double>> global_H_l = global_H;
    vector<vector<double>> global_C_l = global_C;
    vector<double> global_P_l = global_P;
    for (int i = 0; i < global_data.NodesNumber; i++)
        for (int j = 0; j < global_data.NodesNumber; j++)
        {
            global_C_l[i][j] /= global_data.SimulationStepTime;
            global_P_l[i] += global_C_l[i][j] * InitialTemp[j];
            global_H_l[i][j] += global_C_l[i][j];
        }
    global_PC = global_P_l;
    return solveSystem(global_H_l, global_P_l, global_data);
}

vector<double> SOE::solveSystem(vector<vector<double>> &global_H, vector<double> &global_P, Global_data global_data)
{
    vector<double> solution(global_data.NodesNumber, 0.0);
    const int rows = global_H.size();
    const int cols = global_H[0].size();
    for (int i = 0; i < rows; i++)
        global_H[i].push_back(global_P[i]);

    gaussianElimination(global_H);

    for (int i = 0; i < rows; i++)
        solution[i] = global_H[i][cols];
    return solution;
}

void SOE::gaussianElimination(vector<vector<double>> &global_H)
{
    const int rows = global_H.size();
    const int cols = global_H[0].size();

    for (int i = 0; i < rows - 1; i++)
        for (int k = i + 1; k < rows; k++)
        {
            double factor = global_H[k][i] / global_H[i][i];
            for (int j = i; j < cols; j++)
                global_H[k][j] -= factor * global_H[i][j];
        }

    for (int i = rows - 1; i >= 0; i--)
    {
        global_H[i][cols - 1] /= global_H[i][i];
        for (int k = i - 1; k >= 0; k--)
            global_H[k][cols - 1] -= global_H[k][i] * global_H[i][cols - 1];
    }
}