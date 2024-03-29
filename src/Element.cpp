#include "../include/Element.h"

vector<vector<double>> Element::calculateMatrixH(int n, UniversalElement eu, Node nodes[4], int Conductivity, int SpecificHeat, int Density)
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
    GausTable table(n);

    for (int k = 0; k < n * n; k++)
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
            {
                matrixH[k][i][j] = Conductivity * (tablePCX[k][i] * tablePCXT[j][k] + tablePCY[k][i] * tablePCYT[j][k]) * det[k];
                matrixC[k][i][j] = SpecificHeat * Density * (eu.N_nodes[k][i] * euN_nodesT[j][k]) * det[k];
                wynikH[i][j] += matrixH[k][i][j] * table.weights[k % n] * table.weights[floor(k / n)];
                wynikC[i][j] += matrixC[k][i][j] * table.weights[k % n] * table.weights[floor(k / n)];
            }
    // printMatrix(wynikH);
    // cout << endl;
    C = wynikC;
    return wynikH;
}

vector<vector<double>> Element::calculateMatrixHBC(int n, UniversalElement eu, Node nodes[4], int Alfa, int Tot)
{
    vector<vector<double>> HBC(4, vector<double>(4, 0.0));
    vector<vector<vector<double>>> transposedSurfaces = transposeSurfaces(eu);
    double Pf[4] = {0}, detB[4] = {0};
    int bok[4] = {0};
    GausTable table(n);
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
                    Pf[j] += eu.Surface3D[l].N[i][j] * Tot * table.weights[i] * Alfa * detB[l];
                }

    for (int i = 0; i < 4; i++)
        P[i] = Pf[i];

    return HBC;
}