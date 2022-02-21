#ifndef _SOPLPSOLVER_H_
#define _SOPLPSOLVER_H_

#include "SOPMatrix.hpp"
#include "SOPModel.hpp"
#include "SOPpreprocessing.hpp"

class SOPlpSolve:public SOPSolver{
private:
    size_t p_iterations;
    size_t d_iterations;
    size_t n_iterations;
    size_t b_iterations;
protected:

public:
    SOPlpSolve(){

    }

    bool virtual solve(SOPModel &model, int algorithm);
    bool optimality(const SOPVector input);
    int largest_coeff(const SOPVector input);
    bool primal_simplex(SOPModel &model);
    bool dual_simplex(SOPModel &model);
    bool barrier_method();
    bool netword_simplex();
};

#endif