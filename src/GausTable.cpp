#include "../include/UniversalElement.h"

GausTable::GausTable(int n)
{
    if (n == 2)
    {
        nodes = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
        weights = {1.0, 1.0};
    }
    else if (n == 3)
    {
        nodes = {-sqrt(3.0 / 5.0), 0.0, sqrt(3.0 / 5.0)};
        weights = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    }
    else if (n == 4)
    {
        nodes = {-sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0)), -sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(6.0 / 5.0)), sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(6.0 / 5.0))};
        weights = {(18.0 - sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 + sqrt(30.0)) / 36.0, (18.0 - sqrt(30.0)) / 36.0};
    }
    else
    {
        cerr << "Invalid value. Provide the value {2,3,4}" << endl;
    }
}

double GausTable::Gauss_solution(int n, int m)
{
    GausTable table(n);

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