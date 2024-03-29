#include "../include/UniversalElement.h"

UniversalElement::UniversalElement(int n)
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

void UniversalElement::shape_functions()
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

void UniversalElement::N_data(int n)
{
    GausTable table(n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < 4; j++)
        {
            Surface3D[0].N[i][j] = N_function[j](table.nodes[i], -1.0);
            Surface3D[1].N[i][j] = N_function[j](1.0, table.nodes[i]);
            Surface3D[2].N[i][j] = N_function[j](table.nodes[n - 1 - i], 1.0);
            Surface3D[3].N[i][j] = N_function[j](-1.0, table.nodes[n - 1 - i]);
        }
}

void UniversalElement::N_data_nodes(int n)
{
    GausTable table(n);
    for (int i = 0; i < n * n; i++)
        for (int j = 0; j < 4; j++)
            N_nodes[i][j] = N_function[j](table.nodes[i % n], table.nodes[floor(i / n)]);
}

void UniversalElement::print_Surface3D(int n)
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

void UniversalElement::ksi_eta_data(int n)
{
    GausTable table(n);
    for (int i = 0; i < n * n; i++)
        for (int j = 0; j < 4; j++)
        {
            ksi[i][j] = ksi_function[j](table.nodes[floor(i / n)]);
            eta[i][j] = eta_function[j](table.nodes[i % n]);
        }
}

void UniversalElement::print(int n)
{
    for (int i = 0; i < n * n; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << ksi[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    for (int i = 0; i < n * n; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            cout << eta[i][j] << " ";
        }
        cout << endl;
    }
}