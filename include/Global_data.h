#ifndef GLOBAL_DATA_H
#define GLOBAL_DATA_H

#include <iostream>
#include <fstream>
#include <string>
#include "Grid.h"

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

    void readFromFile(const std::string &filename, Grid &grid);

    void print(Grid &grid);
};

#endif