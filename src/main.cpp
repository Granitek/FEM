#include <sstream>
#include "../include/Global_data.h"
#include "../include/GausTable.h"
#include "../include/UniversalElement.h"
#include "../include/SOE.h"
#include "../include/MatrixFunctions.h"
#include "../include/Element.h"
using namespace std;

void writeVTKFile(const string &fileName, Global_data global_data, Grid grid, vector<double> SolutionC)
{
    ofstream vtkFile(fileName);

    if (!vtkFile.is_open())
    {
        cerr << "Failed to open the file: " << fileName << endl;
        return;
    }

    vtkFile << "# vtk DataFile Version 2.0\n";
    vtkFile << "Unstructured Grid Example\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET UNSTRUCTURED_GRID\n\n";

    vtkFile << "POINTS " << global_data.NodesNumber << " float\n";
    for (int i = 0; i < global_data.NodesNumber; ++i)
    {
        vtkFile << grid.Nodes[i].x << " " << grid.Nodes[i].y << " "
                << "0\n";
    }
    vtkFile << "\n";

    vtkFile << "CELLS " << global_data.ElementsNumber << " " << global_data.ElementsNumber * 5 << "\n";
    for (int i = 0; i < global_data.ElementsNumber; ++i)
    {
        vtkFile << "4 ";
        for (int j = 0; j < 4; ++j)
        {
            vtkFile << grid.Elements[i].ID[j] - 1 << " ";
        }
        vtkFile << "\n";
    }
    vtkFile << "\n";

    vtkFile << "CELL_TYPES " << global_data.ElementsNumber << "\n";
    for (int i = 0; i < global_data.ElementsNumber; ++i)
    {
        vtkFile << "9\n";
    }
    vtkFile << "\n";

    vtkFile << "POINT_DATA " << global_data.NodesNumber << "\n";
    vtkFile << "SCALARS Temp float 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < global_data.NodesNumber; ++i)
    {
        vtkFile << SolutionC[i] << "\n";
    }

    vtkFile.close();
    cout << "VTK file has been successfully created: " << fileName << endl;
}

int main()
{
    // Read from test file
    Global_data globalData;
    Grid grid;
    globalData.readFromFile("Test2.txt", grid);
    // globalData.print(grid);

    // n - number of integration points
    int n = 2;

    UniversalElement eu(n);
    eu.ksi_eta_data(n);
    // eu.print(n);

    Node nodes[4];
    eu.N_data(n);
    eu.N_data_nodes(n);
    SOE soe(globalData.NodesNumber);
    // cout << "Surface3D: " << endl;
    // eu.print_Surface3D(n);

    // Calculate matrix H, HBC and P for each Element
    for (Element &e : grid.Elements)
    {
        for (int i = 0; i < 4; i++)
        {
            nodes[i] = grid.Nodes[e.ID[i] - 1];
        }
        e.H = e.calculateMatrixH(n, eu, nodes, globalData.Conductivity, globalData.SpecificHeat, globalData.Density);
        e.HBC = e.calculateMatrixHBC(n, eu, nodes, globalData.Alfa, globalData.Tot);
    }

    // H, HBC and P for each Element
    //  for (Element &e : grid.Elements)
    //  {
    //      cout << "H:" << endl;
    //      printMatrix(e.H);
    //      cout << endl;
    //      cout << "C:" << endl;
    //      printMatrix(e.C);
    //      cout << endl;
    //  }
    // for (Element &e : grid.Elements)
    // {
    //     cout << "HBC:" << endl;
    //     printMatrix(e.HBC);
    //     cout << endl;
    //     cout << "Vector P: ";
    //     for (int i = 0; i < 4; i++)
    //         cout << e.P[i] << " ";
    //     cout << endl
    //          << endl;
    // }

    // Aggregate
    for (Element &e : grid.Elements)
    {
        soe.global_H = soe.aggregate_H_C(e, globalData);
        soe.global_P = soe.aggregate_P(e, globalData);
    }

    // Global H
    // cout << "Global H: " << endl;
    // printMatrix(soe.global_H);
    // cout << endl << "BC: ";

    // Global P
    // for (int i = 0; i < globalData.NodesNumber; i++)
    //  {
    //      cout << soe.global_P[i] << " ";
    //  }
    // cout << endl;

    // Global C
    // cout << endl<< "Global C: " << endl;
    // printMatrix(soe.global_C);
    // cout << endl;

    int SimulationNumbers = globalData.SimulationTime / globalData.SimulationStepTime;
    vector<double> solutionC;
    vector<double> InitialTemp;
    InitialTemp.resize(globalData.NodesNumber, globalData.InitialTemp);
    solutionC.resize(globalData.NodesNumber, 0.0);
    solutionC = soe.Results(globalData, InitialTemp);

    // BC
    //  cout << "BC: ";
    //  for (int i = 0; i < globalData.NodesNumber; i++)
    //  {
    //      cout << soe.global_PC[i] << " ";
    //  }
    //  cout << endl;

    cout << "Solution of the system of equations" << endl;
    for (int i = 0; i < SimulationNumbers; i++)
    {
        solutionC = soe.Results(globalData, InitialTemp);
        auto min = min_element(solutionC.begin(), solutionC.end());
        auto max = max_element(solutionC.begin(), solutionC.end());
        cout << "Time = " << globalData.SimulationStepTime * (i + 1) << " min_t = " << *min << " max_t = " << *max << endl;

        // Results
        //  for (int i = 0; i < globalData.NodesNumber; i++)
        //  {
        //      cout << solutionC[i] << " ";
        //  }

        // vtk files for Paraview testing
        stringstream ss;
        ss << i;
        string vtkName = "Foo" + ss.str() + ".vtk";
        writeVTKFile(vtkName, globalData, grid, solutionC);
        InitialTemp = solutionC;
    }
    cout << endl;
    return 0;
}
