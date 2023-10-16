#include "SOPModel.hpp"
#include "SOPMatrix.hpp"

SOPVariables SOPModel::create_nonbasic(std::string name, char type, int id, double lb, double ub) {
    SOPVariables new_var = SOPVariables(name, type, id, lb, ub);
    this->nonbasic.push_back(new_var);
    this->num_cols++;
    return new_var;
}

void SOPModel::addMax(SOPObjective _obj) {
    this->c = _obj.get_all_coef();
    for (int i = 0; i < _obj.getSize(); ++i) {
        int tmp = this->c.coeff(i);
        this->c(i) = -tmp;
    }
    _obj.clear();
}

void SOPModel::addMin(SOPObjective _obj) {
    this->c = _obj.get_all_coef();
    _obj.clear();
}

SOPRange SOPModel::addLe(SOPExpression& con, double rhs, std::string con_name) {
    this->A.conservativeResize(this->A.rows() + 1, get_num_nonbasic());
    this->N.conservativeResize(this->N.rows() + 1, get_num_nonbasic());

    for (int i = 0; i < get_num_nonbasic(); ++i) {
        this->A(this->A.rows() - 1, con.get_var(i).get_var_id()) = con.get_coef(i);
        this->N(this->N.rows() - 1, con.get_var(i).get_var_id()) = con.get_coef(i);
    }

    this->b.conservativeResize(this->b.size() + 1);
    this->b(this->b.size() - 1) = rhs;

    // create basic variable
    SOPVariables basic = SOPVariables("x" + std::to_string(this->get_num_cols() + this->get_num_rows() + 1), 'R', this->get_num_cols() + this->get_num_rows(), 0, INT_MAX);
    this->basic.push_back(basic);

    // create range
    SOPRange new_con = SOPRange(con.get_all_coef(), con.get_all_vars(), rhs, 'L', con_name);
    this->cons.push_back(new_con);
    this->num_rows++;
    con.clear();

    return new_con;
}

SOPRange SOPModel::addGe(SOPExpression& con, double rhs, std::string con_name) {
    this->A.conservativeResize(this->A.rows() + 1, get_num_nonbasic());
    this->N.conservativeResize(this->N.rows() + 1, get_num_nonbasic());

    for (int i = 0; i < get_num_nonbasic(); ++i) {
        this->A(this->A.rows() - 1, con.get_var(i).get_var_id()) = -con.get_coef(i);
        this->N(this->N.rows() - 1, con.get_var(i).get_var_id()) = -con.get_coef(i);
    }

    this->b.conservativeResize(this->b.size() + 1);
    this->b(this->b.size() - 1) = rhs;

    // create basic variable
    SOPVariables basic = SOPVariables("x" + std::to_string(this->get_num_cols() + this->get_num_rows() + 1), 'R', this->get_num_cols() + this->get_num_rows(), 0, INT_MAX);
    this->basic.push_back(basic);

    // create range
    SOPRange new_con = SOPRange(con.get_all_coef(), con.get_all_vars(), rhs, 'G', con_name);
    this->cons.push_back(new_con);
    this->num_rows++;
    con.clear();

    return new_con;
}

SOPRange SOPModel::addEq(SOPExpression& con, double rhs, std::string con_name) {
    this->A.conservativeResize(this->A.rows() + 1, get_num_nonbasic());
    this->N.conservativeResize(this->N.rows() + 1, get_num_nonbasic());

    for (int i = 0; i < get_num_nonbasic(); ++i) {
        this->A(this->A.rows() - 1, con.get_var(i).get_var_id()) = con.get_coef(i);
        this->N(this->N.rows() - 1, con.get_var(i).get_var_id()) = con.get_coef(i);
    }

    this->b.conservativeResize(this->b.size() + 1);
    this->b(this->b.size() - 1) = rhs;

    // create basic variable
    // SOPVariables basic = SOPVariables("x"+std::to_string(this->get_num_cols()+this->get_num_rows()), 'R'
    //, this->get_num_cols()+this->get_num_rows(), 0, INT_MAX);

    // create range
    SOPRange new_con = SOPRange(con.get_all_coef(), con.get_all_vars(), rhs, 'E', con_name);
    this->cons.push_back(new_con);
    this->num_rows++;
    con.clear();

    return new_con;
}

void SOPModel::create_basis() {
    this->B = SOPMatrix::Identity(this->num_rows, this->num_rows);
    this->A.conservativeResize(this->A.rows(), this->nonbasic.size() + this->basic.size());
    this->A.leftCols(this->N.cols()) = this->N;
    this->A.rightCols(this->B.cols()) = this->B;
}

int SOPModel::get_num_rows() {
    return this->num_rows;
}

int SOPModel::get_num_cols() {
    return this->num_cols;
}

int SOPModel::get_num_nonbasic() {
    return this->nonbasic.size();
}

int SOPModel::get_num_basic() {
    return this->basic.size();
}

int SOPModel::get_all_variables() {
    return this->nonbasic.size() + this->basic.size();
}

SOPVariables SOPModel::get_nonbasic(int position) {
    return this->nonbasic[position];
}

SOPVariables SOPModel::get_basic(int position) {
    return this->basic[position];
}

SOPVector SOPModel::get_header() {
    return this->c;
}

SOPVector SOPModel::get_rhs() {
    return this->b;
}

SOPMatrix SOPModel::get_all_range() {
    return this->A;
}

SOPMatrix SOPModel::get_basis() {
    return this->B;
}

SOPMatrix SOPModel::get_N() {
    return this->N;
}

void SOPModel::set_header(int position, double value) {
    this->c(position) = value;
}

void SOPModel::set_rhs(int position, double value) {
    this->b(position) = value;
}

void SOPModel::set_nonbasic(int position, SOPVariables var) {
    this->nonbasic[position] = var;
}

void SOPModel::set_basic(int position, SOPVariables var) {
    this->basic[position] = var;
}

void SOPModel::set_eta(int leaving_index, std::vector<double> delta_x_B) {
    this->eta.resize(this->num_rows, this->num_rows);
    this->eta.setZero();
    this->eta.diagonal().array() += 1;
    for (int i = 0; i < this->num_rows; ++i) {
        double value = 0;
        if (i != leaving_index)
            value = -(delta_x_B[i] / delta_x_B[leaving_index]);
        else
            value = 1 / delta_x_B[leaving_index];
        this->eta(i, leaving_index) = value;
    }
    std::cout << "create identical matrix" << std::endl;
}

void SOPModel::update_basis() {
    this->B = this->eta * this->B;
}

void SOPModel::clear() {
}

bool SOPSolver::solve(SOPModel model, int algorithm) {
    return true;
}