#ifndef GRID_H
#define GRID_H

#include <vector>
#include <fstream>
#include "Element.h"

using namespace std;

struct Grid
{
    vector<Element> Elements;
    vector<Node> Nodes;

    void readNodesFromFile(ifstream &file, int nodesNumber);
    void readBCFromFile(ifstream &file, Grid &grid);
    void readElementsFromFile(ifstream &file, int elementsNumber);
};

#endif