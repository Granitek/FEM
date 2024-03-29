#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#include <sstream>

using namespace std;

struct Global_data
{
    int SimulationTime;
    int SimulationStepTime;
    int Conductivity;
    int Alfa;
    int Tot;
    int InitialTemp;
    int Density;
    int SpecificHeat;
    int NodesNumber;
    int ElementsNumber;
};

struct Node
{
    double x;
    double y;
    int BC = 0;
};

struct Element
{
    int ID[4];
    vector<vector<double>> H;
    vector<vector<double>> HBC;
    vector<vector<double>> C;
    double P[4];
};

struct Grid
{
    vector<Element> Elements;
    vector<Node> Nodes;
};

struct GausTable
{
    vector<double> nodes;
    vector<double> weights;
};

struct Surface
{
    vector<vector<double>> N;
    Surface(int n)
    {
        N = vector<vector<double>>(n, vector<double>(4));
    }
};

GausTable Gaus_table_n(int n);
double Gauss_solution(int n, int m);

struct UniversalElement
{
    vector<vector<double>> ksi;
    vector<vector<double>> eta;
    vector<Surface> Surface3D;
    vector<vector<double>> N_nodes;

    UniversalElement(int n)
    {
        Surface surface(n);
        Surface3D = vector<Surface>(4, surface);
        if (n * n != 4 && n * n != 9 && n * n != 16)
            throw invalid_argument("Invalid value. Provide the value {2,3,4}");
        ksi = vector<vector<double>>(n * n, vector<double>(4));
        eta = vector<vector<double>>(n * n, vector<double>(4));
        N_nodes = vector<vector<double>>(n * n, vector<double>(4));
        shape_functions();
    }

    vector<function<double(double)>> ksi_function = vector<function<double(double)>>();
    vector<function<double(double)>> eta_function = vector<function<double(double)>>();
    vector<function<double(double, double)>> N_function = vector<function<double(double, double)>>();
    void shape_functions()
    {
        ksi_function.push_back([](double x)
                               { return -0.25 * (1 - x); });
        ksi_function.push_back([](double x)
                               { return 0.25 * (1 - x); });
        ksi_function.push_back([](double x)
                               { return 0.25 * (1 + x); });
        ksi_function.push_back([](double x)
                               { return -0.25 * (1 + x); });

        eta_function.push_back([](double x)
                               { return -0.25 * (1 - x); });
        eta_function.push_back([](double x)
                               { return -0.25 * (1 + x); });
        eta_function.push_back([](double x)
                               { return 0.25 * (1 + x); });
        eta_function.push_back([](double x)
                               { return 0.25 * (1 - x); });

        N_function.push_back([](double ksi, double eta)
                             { return 0.25 * (1 - ksi) * (1 - eta); });
        N_function.push_back([](double ksi, double eta)
                             { return 0.25 * (1 + ksi) * (1 - eta); });
        N_function.push_back([](double ksi, double eta)
                             { return 0.25 * (1 + ksi) * (1 + eta); });
        N_function.push_back([](double ksi, double eta)
                             { return 0.25 * (1 - ksi) * (1 + eta); });
    }

    void N_data(int n)
    {
        GausTable table = Gaus_table_n(n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < 4; j++)
            {
                Surface3D[0].N[i][j] = N_function[j](table.nodes[i], -1.0);
                Surface3D[1].N[i][j] = N_function[j](1.0, table.nodes[i]);
                Surface3D[2].N[i][j] = N_function[j](table.nodes[n - 1 - i], 1.0);
                Surface3D[3].N[i][j] = N_function[j](-1.0, table.nodes[n - 1 - i]);
            }
    }

    void N_data_nodes(int n)
    {
        GausTable table = Gaus_table_n(n);
        for (int i = 0; i < n * n; i++)
            for (int j = 0; j < 4; j++)
                N_nodes[i][j] = N_function[j](table.nodes[i % n], table.nodes[floor(i / n)]);
    }

    void print_Surface3D(int n)
    {
        for (int k = 0; k < 4; k++)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < 4; j++)
                    cout << Surface3D[k].N[i][j] << " ";
                cout << endl;
            }
            cout << endl;
        }
    }

    void ksi_eta_data(int n)
    {
        GausTable table = Gaus_table_n(n);
        for (int i = 0; i < n * n; i++)
            for (int j = 0; j < 4; j++)
            {
                ksi[i][j] = ksi_function[j](table.nodes[floor(i / n)]);
                eta[i][j] = eta_function[j](table.nodes[i % n]);
            }
    }
};

void gaussianElimination(vector<vector<double>> &global_H);
vector<double> solveSystem(vector<vector<double>> &global_H, vector<double> &global_P, Global_data global_data);
void printMatrix(const vector<vector<double>> &matrix);

struct SOE
{
    vector<vector<double>> global_H;
    vector<double> global_P;
    vector<vector<double>> global_C;
    vector<double> global_PC;

    SOE(int NodesNumber)
    {
        global_H.resize(NodesNumber, vector<double>(NodesNumber, 0.0));
        global_P.resize(NodesNumber, 0.0);
        global_C.resize(NodesNumber, vector<double>(NodesNumber, 0.0));
        global_PC.resize(NodesNumber, 0.0);
    }

    vector<vector<double>> agregate_H_C(Element element, Global_data global_data)
    {
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
            {
                global_H[element.ID[i] - 1][element.ID[j] - 1] += element.H[i][j] + element.HBC[i][j];
                global_C[element.ID[i] - 1][element.ID[j] - 1] += element.C[i][j];
            }
        return global_H;
    }

    vector<double> agregate_P(Element element, Global_data global_data)
    {
        for (int i = 0; i < 4; i++)
            global_P[element.ID[i] - 1] += element.P[i];
        return global_P;
    }

    vector<double> Results(Global_data global_data, vector<double> InitialTemp)
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
};

void readNodesFromFile(ifstream &file, vector<Node> &nodes, int nodesNumber);
void readElementsFromFile(ifstream &file, vector<Element> &elements, int elementsNumber);
void readGlobalDataFromFile(const string &filename, Global_data &globalData, Grid &grid);

vector<vector<double>> transposeMatrix(const vector<vector<double>> &matrix);

vector<vector<double>> calculateMatrixH(int n, UniversalElement eu, Node nodes[4], Global_data global_data, Element &e);

void readBCFromFile(ifstream &file, Grid &grid);
vector<vector<vector<double>>> transposeSurfaces(const UniversalElement &eu);
vector<vector<double>> calculateHBC(int n, UniversalElement eu, Node nodes[4], int Alfa, int Tot, Element &e);

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
    cout << "----------------------------------LAB 1----------------------------------" << endl;
    Global_data globalData;
    Grid grid;
    readGlobalDataFromFile("Test2.txt", globalData, grid);
    int n = 2;

    // cout << "SimulationTime: " << globalData.SimulationTime << endl;
    // cout << "SimulationStepTime: " << globalData.SimulationStepTime << endl;
    // cout << "Conductivity: " << globalData.Conductivity << endl;
    // cout << "Alfa: " << globalData.Alfa << endl;
    // cout << "Tot: " << globalData.Tot << endl;
    // cout << "InitialTemp: " << globalData.InitialTemp << endl;
    // cout << "Density: " << globalData.Density << endl;
    // cout << "SpecificHeat: " << globalData.SpecificHeat << endl;
    // cout << "NodesNumber: " << globalData.NodesNumber << endl;
    // cout << "ElementsNumber: " << globalData.ElementsNumber << endl;

    // cout << "Nodes:" << endl;
    // for (const auto &node : grid.Nodes)
    //     cout << node.x << ", " << node.y << " "
    //          << node.BC << endl;

    // cout << "Elements:" << endl;
    // for (const auto &element : grid.Elements)
    //     cout << element.ID[0] << ", " << element.ID[1] << ", " << element.ID[2] << ", " << element.ID[3] << endl;
    cout << "----------------------------------LAB 2----------------------------------" << endl;
    // cout << "Dla n = 2: " << Gauss_solution(2, 1) << endl
    //      << Gauss_solution(3, 1) << endl
    //      << Gauss_solution(3, 2) << endl
    //      << Gauss_solution(3, 2) << endl;
    cout << "----------------------------------LAB 3----------------------------------" << endl;
    UniversalElement eu(n);
    eu.ksi_eta_data(n);
    // for (int i = 0; i < n * n; i++)
    // {
    //     for (int j = 0; j < 4; j++)
    //     {
    //         cout << eu.ksi[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << endl;
    // for (int i = 0; i < n * n; i++)
    // {
    //     for (int j = 0; j < 4; j++)
    //     {
    //         cout << eu.eta[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    cout << "----------------------------------LAB 4----------------------------------" << endl;
    Node nodes[4];
    eu.N_data(n);
    eu.N_data_nodes(n);
    SOE soe(globalData.NodesNumber);
    // cout << "Surface3D: " << endl;
    // eu.print_Surface3D(n);
    for (Element &e : grid.Elements)
    {
        for (int i = 0; i < 4; i++)
        {
            nodes[i] = grid.Nodes[e.ID[i] - 1];
        }
        e.H = calculateMatrixH(n, eu, nodes, globalData, e);
        e.HBC = calculateHBC(n, eu, nodes, globalData.Alfa, globalData.Tot, e);
    }
    // for (Element &e : grid.Elements)
    // {
    //     cout << "H:" << endl;
    //     printMatrix(e.H);
    //     cout << endl;
    //     cout << "C:" << endl;
    //     printMatrix(e.C);
    //     cout << endl;
    // }

    // for (Element &e : grid.Elements)
    // {
    //     cout << "HBC:" << endl;
    //     printMatrix(e.HBC);
    //     cout << endl;
    //     cout << "Wektor P: ";
    //     for (int i = 0; i < 4; i++)
    //         cout << e.P[i] << " ";
    //     cout << endl
    //          << endl;
    // }

    for (Element &e : grid.Elements)
    {
        soe.global_H = soe.agregate_H_C(e, globalData);
        soe.global_P = soe.agregate_P(e, globalData);
    }

    cout << "----------------------------------LAB 6----------------------------------" << endl;
    cout << "Global H: " << endl;
    printMatrix(soe.global_H);
    cout << endl
         << "BC: ";
    for (int i = 0; i < globalData.NodesNumber; i++)
    {
        cout << soe.global_P[i] << " ";
    }
    cout << endl;
    // vector<double> solution = solveSystem(soe.global_H, soe.global_P);
    // cout << "Rozwiązanie układu równań:" << endl;
    // for (int i = 0; i < globalData.NodesNumber; i++)
    // {
    //     cout << solution[i] << " ";
    // }
    cout << endl
         << "Global C: " << endl;
    printMatrix(soe.global_C);
    cout << endl;

    int SimulationNumbers = globalData.SimulationTime / globalData.SimulationStepTime;
    vector<double> solutionC;
    vector<double> InitialTemp;
    InitialTemp.resize(globalData.NodesNumber, globalData.InitialTemp);
    solutionC.resize(globalData.NodesNumber, 0.0);
    solutionC = soe.Results(globalData, InitialTemp);

    cout << "BC: ";
    for (int i = 0; i < globalData.NodesNumber; i++)
    {
        cout << soe.global_PC[i] << " ";
    }
    cout << endl;
    cout << "Solution of the system of equations" << endl;
    for (int i = 0; i < SimulationNumbers; i++)
    {
        solutionC = soe.Results(globalData, InitialTemp);
        auto min = min_element(solutionC.begin(), solutionC.end());
        auto max = max_element(solutionC.begin(), solutionC.end());
        cout << "Time = " << globalData.SimulationStepTime * (i + 1) << " min_t = " << *min << " max_t = " << *max << endl;
        // for (int i = 0; i < globalData.NodesNumber; i++)
        // {
        //     cout << solutionC[i] << " ";
        // }
        stringstream ss;
        ss << i;
        string vtkName = "Foo" + ss.str() + ".vtk";
        writeVTKFile(vtkName, globalData, grid, solutionC);
        InitialTemp = solutionC;
    }
    cout << endl;
    return 0;
}

void readGlobalDataFromFile(const string &filename, Global_data &globalData, Grid &grid)
{
    ifstream file(filename);
    if (!file.is_open())
    {
        cerr << "Failed to open the file: " << filename << endl;
        return;
    }
    string key, a;
    while (file >> key)
    {
        if (key == "SimulationTime")
            file >> globalData.SimulationTime;
        else if (key == "SimulationStepTime")
            file >> globalData.SimulationStepTime;
        else if (key == "Conductivity")
            file >> globalData.Conductivity;
        else if (key == "Alfa")
            file >> globalData.Alfa;
        else if (key == "Tot")
            file >> globalData.Tot;
        else if (key == "InitialTemp")
            file >> globalData.InitialTemp;
        else if (key == "Density")
            file >> globalData.Density;
        else if (key == "SpecificHeat")
            file >> globalData.SpecificHeat;
        else if (key == "Nodes")
        {
            file >> a;
            file >> globalData.NodesNumber;
        }
        else if (key == "Elements")
        {
            file >> a;
            file >> globalData.ElementsNumber;
        }
        else if (key == "*Node")
            readNodesFromFile(file, grid.Nodes, globalData.NodesNumber);
        else if (key == "*Element,")
            readElementsFromFile(file, grid.Elements, globalData.ElementsNumber);
        else if (key == "*BC")
            readBCFromFile(file, grid);
    }

    file.close();
}
void readElementsFromFile(ifstream &file, vector<Element> &elements, int elementsNumber)
{
    elements.resize(elementsNumber);
    int id, n1, n2, n3, n4;
    char coma;
    string line;
    file >> line;
    for (int i = 0; i < elementsNumber; ++i)
    {
        file >> id >> coma >> n1 >> coma >> n2 >> coma >> n3 >> coma >> n4;
        elements[i].ID[0] = n1;
        elements[i].ID[1] = n2;
        elements[i].ID[2] = n3;
        elements[i].ID[3] = n4;
    }
}
void readNodesFromFile(ifstream &file, vector<Node> &nodes, int nodesNumber)
{
    nodes.resize(nodesNumber);
    int id;
    char coma;
    double x, y;
    for (int i = 0; i < nodesNumber; ++i)
    {
        file >> id >> coma >> x >> coma >> y;
        nodes[i] = {x, y};
    }
}

GausTable Gaus_table_n(int n)
{
    GausTable table;

    if (n == 2)
    {
        table.nodes = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
        table.weights = {1.0, 1.0};
    }
    else if (n == 3)
    {
        table.nodes = {-sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0)};
        table.weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    }
    else if (n == 4)
    {
        table.nodes = {-sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0))};
        table.weights = {(18.0 - sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0};
    }
    else
    {
        cerr << "Invalid value. Provide the value {2,3,4}" << endl;
    }

    return table;
}
double Gauss_solution(int n, int m)
{
    GausTable table = Gaus_table_n(n);

    double integral, x, fx, y;
    if (m == 1)
    {
        for (int i = 0; i < n; ++i)
        {
            x = table.nodes[i];
            // cout << x << endl;
            fx = 5.0 * x * x + 3.0 * x + 6.0;
            integral += table.weights[i] * fx;
        }
    }
    else if (m == 2)
    {
        for (int i = 0; i < n; i++)
        {
            x = table.nodes[i];
            for (int j = 0; j < n; j++)
            {
                y = table.nodes[j];
                fx = 5.0 * x * x * y * y + 3.0 * x * y + 6.0;
                integral += table.weights[i] * table.weights[j] * fx;
            }
        }
    }
    return integral;
}

vector<vector<double>> transposeMatrix(const vector<vector<double>> &matrix)
{
    int rows = matrix.size();
    int cols = matrix[0].size();

    vector<vector<double>> transposedArray(cols, vector<double>(rows, 0.0));

    for (int i = 0; i < cols; ++i)
    {
        for (int j = 0; j < rows; ++j)
        {
            transposedArray[i][j] = matrix[j][i];
        }
    }

    return transposedArray;
}
void printMatrix(const vector<vector<double>> &matrix)
{
    for (const auto &row : matrix)
    {
        for (double element : row)
        {
            cout << element << " ";
        }
        cout << endl;
    }
}
vector<vector<double>> calculateMatrixH(int n, UniversalElement eu, Node nodes[4], Global_data global_data, Element &e)
{
    vector<double> dxdksi(n * n, 0.0);
    vector<double> dydksi(n * n, 0.0);
    vector<double> dxdeta(n * n, 0.0);
    vector<double> dydeta(n * n, 0.0);
    for (int i = 0; i < n * n; i++)
        for (int j = 0; j < 4; j++)
        {
            dxdksi[i] += eu.ksi[i][j] * nodes[j].x;
            dydksi[i] += eu.ksi[i][j] * nodes[j].y;
            dxdeta[i] += eu.eta[i][j] * nodes[j].x;
            dydeta[i] += eu.eta[i][j] * nodes[j].y;
        }
    vector<double> det(n * n, 0.0);
    for (int i = 0; i < n * n; i++)
    {
        // LICZENIE JAKOBIANU
        det[i] = dxdksi[i] * dydeta[i] - dydksi[i] * dxdeta[i];
        det[i] = 1 / det[i];
        // ODWRÓCENIE JAKOBIANU
        dxdksi[i] *= det[i];
        dydksi[i] *= -1 * det[i];
        dxdeta[i] *= -1 * det[i];
        dydeta[i] *= det[i];

        // ODWRÓCENIE MACIERZY
        double pom;
        pom = dxdksi[i];
        dxdksi[i] = dydeta[i];
        dydeta[i] = pom;
        pom = dydksi[i];
        dydksi[i] = dxdeta[i];
        dxdeta[i] = pom;

        det[i] = 1 / det[i];
    }

    // for (int i = 0; i < n * n; i++)
    // {
    //     cout << "Dla pc" << i << ":" << endl
    //          << dxdksi[i] << " " << dydksi[i] << endl
    //          << dxdeta[i] << " " << dydeta[i] << endl;
    // }
    vector<vector<double>> tablePCX(n * n, vector<double>(4, 0.0));
    vector<vector<double>> tablePCY(n * n, vector<double>(4, 0.0));
    // cout << endl
    //      << "tablePCX:" << endl;
    for (int i = 0; i < n * n; i++)
        for (int j = 0; j < 4; j++)
        {
            tablePCX[i][j] = dxdksi[i] * eu.ksi[i][j] + dxdeta[i] * eu.eta[i][j];
            tablePCY[i][j] = dydksi[i] * eu.ksi[i][j] + dydeta[i] * eu.eta[i][j];
        }

    // printMatrix(tablePCX);
    // cout << endl
    //      << "tablePCY:" << endl;
    // printMatrix(tablePCY);
    // cout << endl;

    vector<vector<vector<double>>> matrixH(n * n, vector<vector<double>>(4, vector<double>(4, 0.0)));
    vector<vector<vector<double>>> matrixC(n * n, vector<vector<double>>(4, vector<double>(4, 0.0)));
    vector<vector<double>> tablePCXT = transposeMatrix(tablePCX);
    vector<vector<double>> tablePCYT = transposeMatrix(tablePCY);
    vector<vector<double>> euN_nodesT = transposeMatrix(eu.N_nodes);
    vector<vector<double>> wynikH(4, vector<double>(4, 0.0));
    vector<vector<double>> wynikC(4, vector<double>(4, 0.0));
    GausTable table = Gaus_table_n(n);

    for (int k = 0; k < n * n; k++)
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
            {
                matrixH[k][i][j] = global_data.Conductivity * (tablePCX[k][i] * tablePCXT[j][k] + tablePCY[k][i] * tablePCYT[j][k]) * det[k];
                matrixC[k][i][j] = global_data.SpecificHeat * global_data.Density * (eu.N_nodes[k][i] * euN_nodesT[j][k]) * det[k];
                wynikH[i][j] += matrixH[k][i][j] * table.weights[k % n] * table.weights[floor(k / n)];
                wynikC[i][j] += matrixC[k][i][j] * table.weights[k % n] * table.weights[floor(k / n)];
            }
    // printMatrix(wynikH);
    // cout << endl;
    e.C = wynikC;
    return wynikH;
}

void readBCFromFile(ifstream &file, Grid &grid)
{
    int elementID;
    string coma;
    while (file >> elementID)
    {
        file >> coma;
        grid.Nodes[elementID - 1].BC = 1;
    }
}
vector<vector<vector<double>>> transposeSurfaces(const UniversalElement &eu)
{
    vector<vector<vector<double>>> transposedSurfaces;
    for (const Surface &surface : eu.Surface3D)
    {
        transposedSurfaces.push_back(transposeMatrix(surface.N));
    }
    return transposedSurfaces;
}
vector<vector<double>> calculateHBC(int n, UniversalElement eu, Node nodes[4], int Alfa, int Tot, Element &e)
{
    vector<vector<double>> HBC(4, vector<double>(4, 0.0));
    vector<vector<vector<double>>> transposedSurfaces = transposeSurfaces(eu);
    double P[4] = {0}, detB[4] = {0};
    int bok[4] = {0};
    GausTable table = Gaus_table_n(n);
    for (int i = 0; i < 4; i++)
    {
        if (nodes[i % 4].BC != 0 && nodes[(i + 1) % 4].BC == nodes[i % 4].BC)
            bok[i] = 1;
        detB[i % 4] = sqrt(pow(nodes[(i + 1) % 4].x - nodes[i % 4].x, 2) + pow(nodes[(i + 1) % 4].y - nodes[i % 4].y, 2)) / 2;
    }

    // cout << "Surface3DT: " << endl;
    // for (const auto &matrix : transposedSurfaces)
    // {
    //     for (const auto &row : matrix)
    //     {
    //         for (double val : row)
    //             cout << val << " ";
    //         cout << endl;
    //     }
    //     cout << endl;
    // }

    for (int l = 0; l < 4; l++)
        for (int i = 0; i < n; i++)
            for (int j = 0; j < 4; j++)
                if (bok[l] == 1)
                {
                    for (int k = 0; k < 4; k++)
                    {
                        HBC[j][k] += eu.Surface3D[l].N[i][j] * transposedSurfaces[l][k][i] * Alfa * detB[l] * table.weights[i];
                    }
                    P[j] += eu.Surface3D[l].N[i][j] * Tot * table.weights[i] * Alfa * detB[l];
                }

    for (int i = 0; i < 4; i++)
        e.P[i] = P[i];

    return HBC;
}

void gaussianElimination(vector<vector<double>> &global_H)
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
vector<double> solveSystem(vector<vector<double>> &global_H, vector<double> &global_P, Global_data global_data)
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