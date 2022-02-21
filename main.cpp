#include "SOPMatrix.hpp"
#include "SOPModel.hpp"
#include "SOPlpsolver.hpp"
#include "SOPqpsolver.hpp"
#include "SOPpreprocessing.hpp"
#include "SOPbcsolver.hpp"

int main(){
    std::cout<<"hello, world"<<std::endl;
    #pragma region testing for Eigen Lib
    SOPMatrix m1;
    m1 = SOPMatrix::Constant(2,3,1);
    std::cout<<m1<<std::endl;
    #pragma endregion

    #pragma region Mathematical Programming Modeling
    SOPModel model;
    SOPExpression expr;

    // create decision variables
    SOPVariables x1 = model.create_nonbasic("x1", 'F', 0, 0, INT_MAX);
    SOPVariables x2 = model.create_nonbasic("x2", 'F', 1, 0, INT_MAX);
    SOPVariables x3 = model.create_nonbasic("x3", 'F', 2, 0, INT_MAX);

    std::cout<<x1.get_var_name()<<std::endl;
    std::cout<<x2.get_var_name()<<std::endl;

    // create objective function
    SOPObjective obj;
    obj.addTerm(3, x1);
    obj.addTerm(6, x2);
    obj.addTerm(0, x3);
    model.addMin(obj);

    std::cout<<model.get_header()<<std::endl;

    // create bunch of constraints
    expr.addTerm(-1, x1);
    expr.addTerm(2, x2);
    expr.addTerm(-1, x3);
    model.addLe(expr, -6, "cons1");

    expr.addTerm(-2, x1);
    expr.addTerm(-1, x2);
    expr.addTerm(2, x3);
    model.addLe(expr, -8, "cons2");

    /*expr.addTerm(0, x1);
    expr.addTerm(1, x2);
    model.addLe(expr, 5, "cons2");*/

    model.create_basis();
    std::cout<<model.get_basis()<<std::endl;
    std::cout<<model.get_all_range()<<std::endl;
    std::cout<<model.get_N()<<std::endl;
    std::cout<<model.get_rhs()<<std::endl;

    for (int i = 0; i < 2; i++) std::cout<<model.get_basic(i).get_var_id()<<std::endl;

    #pragma endregion

    #pragma region solve linear program
    SOPlpSolve lpsolver;
    lpsolver.solve(model, 2);
    std::cout<<model.get_rhs()<<std::endl;
    #pragma endregion

    std::cout<<"finish"<<std::endl;
    return 0;
}
