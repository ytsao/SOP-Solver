# SOP-Solver
The own mathematical programming solver, Search Optimal Points Solver.

I would like to develop my own mathematical programming solver.
For implementing issue, I used C++ with Eigen Library to do complicate linear algebra more efficient.

The ultimate goal is to build a geneal MIP Solver.

In current stage, only include "primal simplex method" & "dual simplex method"
After finishing "barrier method" & "network simlex method", I will move on to implement "The Wolfes' simplex method" for solving quadratic programs.

During the developing, I still study some classical textbook to strengthen my professional knowledge about discrete & combinatorial opitmization.

## Dependency
1. Eigen3+
```cmd
git clone https://gitlab.com/libeigen/eigen.git
cd eigen
mkdir build
cd build
cmake ..
sudo make install
sudo cp -r /usr/local/include/eigen3/Eigen usr/local/include
```	

## How to use?
```cmd
./run.sh
```