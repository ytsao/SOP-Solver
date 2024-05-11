#! /bin/bash
g++ main.cpp ./include/SOPMatrix.cpp ./include/SOPModel.cpp ./include/SOPlpsolver.cpp -o main
./main
