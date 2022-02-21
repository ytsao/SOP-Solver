#ifndef _SOPMODEL_H_
#define _SOPMODEL_H_

#include "SOPMatrix.hpp"

class SOPModel{
private:
    int num_rows=0;
    int num_cols=0;
    std::vector<SOPVariables> nonbasic;
    std::vector<SOPVariables> basic;
    SOPObjective obj;
    std::vector<SOPRange> cons;
    SOPVector b;
    SOPVector c;
    SOPMatrix A;
    SOPMatrix N;
    SOPMatrix B;
    SOPMatrix eta; 
protected:

public:
    SOPModel(){

    }
    
    SOPVariables create_nonbasic(std::string name, char type, int id, double lb, double ub);
    void addMax(SOPObjective _obj);
    void addMin(SOPObjective _obj);
    SOPRange addLe(SOPExpression &con, double rhs, std::string con_name);
    SOPRange addGe(SOPExpression &con, double rhs, std::string con_name);
    SOPRange addEq(SOPExpression &con, double rhs, std::string con_name);
    void create_basis();
    int get_num_rows();
    int get_num_cols();
    int get_num_nonbasic();
    int get_num_basic();
    int get_all_variables();
    SOPVariables get_nonbasic(int position);
    SOPVariables get_basic(int position);
    SOPVector get_header();
    SOPVector get_rhs();
    SOPMatrix get_all_range();
    SOPMatrix get_basis();
    SOPMatrix get_N();
    void set_header(int position, double value);
    void set_rhs(int position, double value);
    void set_nonbasic(int position, SOPVariables var);
    void set_basic(int position, SOPVariables var);
    void set_eta(int leaving_index, std::vector<double> delta_x_B);
    void update_basis();
    void clear();
};

class SOPSolver{
private:

protected:

public:
    SOPSolver(){

    }

    bool solve(SOPModel model, int algorithm);
};

#endif