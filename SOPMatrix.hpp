#ifndef _SOPMATRIX_H_
#define _SOPMATRIX_H_

#include<bits/stdc++.h>
#include<Eigen/Dense>

typedef Eigen::VectorXd SOPVector;
typedef Eigen::MatrixXd SOPMatrix;

class SOPVariables{
private:
    std::string name;
    char type;
    int id;
    double lb;
    double ub;
protected:

public:
    SOPVariables(){

    }

    SOPVariables(std::string _name, char _type, int _id);
    SOPVariables(std::string _name, char _type, int _id, double _lb, double _ub);
    std::string get_var_name();
    int get_var_id();
};

class SOPExpression{
private:
    SOPVector coef;
    std::vector<SOPVariables> vars;
protected:

public:
    SOPExpression(){
        
    }
   
    void addTerm(const double coef, const SOPVariables var);
    void addTerms(const double* coefs, const SOPVariables* vars);
    int getSize();
    double get_coef(const int position);
    SOPVector get_all_coef();
    SOPVariables get_var(const int position);
    std::vector<SOPVariables> get_all_vars();
    void clear();
};

class SOPObjective:public SOPExpression{
private:
    SOPVector coef;
    std::vector<SOPVariables> vars;
protected:

public:
    SOPObjective(){

    }
};

class SOPRange:public SOPExpression{
private:
    SOPVector coef;
    std::vector<SOPVariables> vars;
    double rhs;
    char direction;
    std::string con_name;
protected:

public:
    SOPRange(){

    }

    SOPRange(SOPVector _coef, std::vector<SOPVariables> _vars, double _rhs, char _direction, std::string _con_name);
};

#endif