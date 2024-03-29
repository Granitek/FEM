#include "../include/Global_data.h"
#include "../include/Grid.h"

void Global_data::readFromFile(const std::string &filename, Grid &grid)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Failed to open the file: " << filename << std::endl;
        return;
    }
    std::string key, a;
    while (file >> key)
    {
        if (key == "SimulationTime")
            file >> SimulationTime;
        else if (key == "SimulationStepTime")
            file >> SimulationStepTime;
        else if (key == "Conductivity")
            file >> Conductivity;
        else if (key == "Alfa")
            file >> Alfa;
        else if (key == "Tot")
            file >> Tot;
        else if (key == "InitialTemp")
            file >> InitialTemp;
        else if (key == "Density")
            file >> Density;
        else if (key == "SpecificHeat")
            file >> SpecificHeat;
        else if (key == "Nodes")
        {
            file >> a;
            file >> NodesNumber;
        }
        else if (key == "Elements")
        {
            file >> a;
            file >> ElementsNumber;
        }
        else if (key == "*Node")
            grid.readNodesFromFile(file, NodesNumber);
        else if (key == "*Element,")
            grid.readElementsFromFile(file, ElementsNumber);
        else if (key == "*BC")
            grid.readBCFromFile(file, grid);
    }

    file.close();
}

void Global_data::print(Grid &grid)
{
    cout << "SimulationTime: " << SimulationTime << endl;
    cout << "SimulationStepTime: " << SimulationStepTime << endl;
    cout << "Conductivity: " << Conductivity << endl;
    cout << "Alfa: " << Alfa << endl;
    cout << "Tot: " << Tot << endl;
    cout << "InitialTemp: " << InitialTemp << endl;
    cout << "Density: " << Density << endl;
    cout << "SpecificHeat: " << SpecificHeat << endl;
    cout << "NodesNumber: " << NodesNumber << endl;
    cout << "ElementsNumber: " << ElementsNumber << endl;

    cout << "Nodes:" << endl;
    for (const auto &node : grid.Nodes)
        cout << node.x << ", " << node.y << " "
             << node.BC << endl;

    cout << "Elements:" << endl;
    for (const auto &element : grid.Elements)
        cout << element.ID[0] << ", " << element.ID[1] << ", " << element.ID[2] << ", " << element.ID[3] << endl;
}