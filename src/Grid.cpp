#include "../include/Grid.h"

void Grid::readNodesFromFile(ifstream &file, int nodesNumber)
{
    Nodes.resize(nodesNumber);
    int id;
    char coma;
    double x, y;
    for (int i = 0; i < nodesNumber; ++i)
    {
        file >> id >> coma >> x >> coma >> y;
        Nodes[i] = {x, y};
    }
}

void Grid::readBCFromFile(ifstream &file, Grid &grid)
{
    int elementID;
    string coma;
    while (file >> elementID)
    {
        file >> coma;
        grid.Nodes[elementID - 1].BC = 1;
    }
}

void Grid::readElementsFromFile(ifstream &file, int elementsNumber)
{
    Elements.resize(elementsNumber);
    int id, n1, n2, n3, n4;
    char coma;
    string line;
    file >> line;
    for (int i = 0; i < elementsNumber; ++i)
    {
        file >> id >> coma >> n1 >> coma >> n2 >> coma >> n3 >> coma >> n4;
        Elements[i].ID[0] = n1;
        Elements[i].ID[1] = n2;
        Elements[i].ID[2] = n3;
        Elements[i].ID[3] = n4;
    }
}