## What is it for?
This project simulates heat transfer in a given element. We have the capability to analyze temperature distributions and heat fluxes throughout the structure.<br />
The program produces VTK files, which can be utilized for result analysis in Paraview.

## Setup
Navigate to src folder:
```
$ cd src
```

To run, use the commands:
```
$ make
$ ./main
```
If you want to change the number of integration points or the text file from which the program reads data, modify variables in the main.cpp file.

```
globalData.readFromFile("Test2.txt", grid);
int n = 2;
```
In this case program will read data from Test2.txt and use 2 integration points.