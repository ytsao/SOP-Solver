#include "SOPMatrix.hpp"

SOPVariables::SOPVariables(std::string _name, char _type, int _id) {
    this->name = _name;
    this->type = _type;
    this->id = _id;
    this->lb = 0;
    this->ub = INT_MAX;
}

SOPVariables::SOPVariables(std::string _name, char _type, int _id, double _lb, double _ub) {
    this->name = _name;
    this->type = _type;
    this->id = _id;
    this->lb = lb;
    this->ub = ub;
}

std::string SOPVariables::get_var_name() {
    return this->name;
}

int SOPVariables::get_var_id() {
    return this->id;
}

void SOPExpression::addTerm(const double coef, const SOPVariables var) {
    /*SOPVector tmp = Eigen::Map<SOPVector>(this->coef.data(), this->coef.size());
    SOPVector new_value(1);
    new_value.fill(coef);*/
    // another way to add new value into SOPVector

    this->coef.conservativeResize(this->coef.size() + 1);
    this->coef(this->coef.size() - 1) = coef;

    this->vars.push_back(var);
}

void SOPExpression::addTerms(const double* coef, const SOPVariables* vars) {
    int size = sizeof(coef) / sizeof(double);
    for (int i = 0; i < size; ++i) {
        this->coef.conservativeResize(this->coef.size() + 1);
        this->coef(this->coef.size() - 1) = coef[i];

        this->vars.push_back(vars[i]);
    }
}

int SOPExpression::getSize() {
    return this->vars.size();
}

double SOPExpression::get_coef(const int position) {
    return this->coef.coeff(position);
}

SOPVector SOPExpression::get_all_coef() {
    return this->coef;
}

SOPVariables SOPExpression::get_var(const int position) {
    return this->vars[position];
}

std::vector<SOPVariables> SOPExpression::get_all_vars() {
    return this->vars;
}

void SOPExpression::clear() {
    this->coef.resize(0);
    this->vars.clear();
    this->vars.shrink_to_fit();
}

SOPRange::SOPRange(SOPVector _coef, std::vector<SOPVariables> _vars, double _rhs, char _direction, std::string _con_name) {
    this->coef = _coef;
    this->vars = _vars;
    this->rhs = _rhs;
    this->direction = _direction;
    this->con_name = _con_name;
}