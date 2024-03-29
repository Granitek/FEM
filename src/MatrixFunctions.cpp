#include "../include/MatrixFunctions.h"

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

vector<vector<vector<double>>> transposeSurfaces(const UniversalElement &eu)
{
    vector<vector<vector<double>>> transposedSurfaces;
    for (const Surface &surface : eu.Surface3D)
    {
        transposedSurfaces.push_back(transposeMatrix(surface.N));
    }
    return transposedSurfaces;
}