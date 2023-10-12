#include "SOPlpsolver.hpp"

bool SOPlpSolve::optimality(const SOPVector input){
    for(int i = 0; i < input.size(); ++i){
        if (input.coeff(i) <= 0)return false;
    }
    return true;
}

int SOPlpSolve::largest_coeff(const SOPVector input){
    double benchmark;
    int index = -1;
    benchmark = INT_MAX;
    for(int i = 0; i < input.size(); ++i){
        if (input.coeff(i) < benchmark){
            benchmark = input.coeff(i);
            index = i;
        }
    }
    return index;
}

bool SOPlpSolve::primal_simplex(SOPModel &model){
    // step1. optimality test
    while(!optimality(model.get_header())){
        std::cout<<"finish step 1"<<std::endl;
        
        // step 2. select entering variable
        int entering_index = largest_coeff(model.get_header());
        SOPVariables entering_variable = model.get_nonbasic(entering_index);
        std::cout<<"entering index = "<<entering_index<<std::endl;
        std::cout<<"finish step 2"<<std::endl;

        // step 3. compute primal step direction
        // delta_x_B = B_inverse * N * e_j
        // N * e_j -> take jth column from N_matrix
        // take whole B_matrix
        std::vector<double> delta_x_B;
        for(int i = 0;i < model.get_num_rows(); ++i){
            double value = 0;
            for (int j = 0; j < model.get_num_rows(); ++j)
                value += model.get_basis()(i,j) * model.get_N()(j,entering_index);
            delta_x_B.push_back(value);
        }
        std::cout<<"finish step 3"<<std::endl;

        // step 4. compute primal step length
        // minimum ratio test
        double t = INT_MIN;
        int leaving_index = -1;
        for (int i = 0; i < model.get_num_rows(); ++i){
            double ratio = delta_x_B[i] / model.get_rhs().coeff(i);
            if (ratio > t && ratio != 0){
                t = ratio;
                leaving_index = i;
            }
        }
        std::cout<<"finish step 4"<<std::endl;

        // step 5. select leaving variable
        if (leaving_index == -1) return false;
        t = 1/t;
        SOPVariables leaving_variable = model.get_basic(leaving_index);
        std::cout<<"finish step 5"<<std::endl;

        // step 6. compute dual step direction
        // delta_z_N = -(B_inverse * N_matrix)_transpose * e_i
        // take ith row from B_matrix and whole N_matrix
        std::vector<double> delta_z_N;
        for (int i = 0; i < model.get_num_cols(); ++i){
            double value = 0;
            for (int j = 0; j < model.get_num_rows(); ++j)
                value += model.get_basis()(leaving_index,j) * model.get_N()(j,i);
            delta_z_N.push_back(-value);
        }
        std::cout<<"finish step 6"<<std::endl;

        // step 7. compute dual step length
        // s = z_star_j / delta_z_j
        double s = model.get_header()(entering_index) / delta_z_N[entering_index];
        std::cout<<"finish step 7"<<std::endl;

        // step 8. update current primal and dual solution
        for (int i = 0; i < model.get_num_rows(); ++i){
            model.set_rhs(i, model.get_rhs()(i)-t*delta_x_B[i]);
            // model.set_header(i, model.get_header()(i)-s*delta_z_N[i]);
        }

        for (int i = 0; i < model.get_num_cols(); ++i){
            model.set_header(i, model.get_header()(i)-s*delta_z_N[i]);
        }
        std::cout<<"finish step 8"<<std::endl;

        // step 9. update basis
        // use eta_matrix to update B_matrix
        // eta_inverse_matrix = I -((delta_x_B - e_i)e_i_T) / delta_x_i
        // new_B_inverse_matrix = eta_inverse_matrix * old_B_inverse_matrix
        model.set_eta(leaving_index, delta_x_B);
        std::cout<<"finish set_eta function"<<std::endl;
        model.update_basis();
        std::cout<<"finish update basis"<<std::endl;

        // update N_matrix
        for (int i = 0; i < model.get_num_rows(); ++i) 
            model.get_N()(i,entering_index) = model.get_all_range()(i,leaving_variable.get_var_id());
        std::cout<<"finish N matrix"<<std::endl;

        // update nonbasic and basic set
        model.set_basic(leaving_index, entering_variable);
        model.set_nonbasic(entering_index, leaving_variable);

        model.set_rhs(leaving_index, t);
        model.set_header(entering_index, s);
        std::cout<< "entering variable index:" << entering_index << ", name: " << entering_variable.get_var_name() << std::endl;
        std::cout<< "leaving variable index:" << leaving_index << ", name: " << leaving_variable.get_var_name() << std::endl;
        std::cout<<"finish step 9"<<std::endl;

        this->p_iterations++;
    }
    return true;
}

bool SOPlpSolve::dual_simplex(SOPModel &model){
    // step 1. optimality test
    while(!optimality(model.get_rhs())){
        std::cout<<"finish step 1"<<std::endl;

        // step 2. select leaving variables
        // pick i in (i in B:x_start_B < 0)
        int leaving_index = largest_coeff(model.get_rhs());
        SOPVariables leaving_variable = model.get_basic(leaving_index);
        std::cout<<"finish step 2"<<std::endl;

        // step 3. compute dual step direction
        // delta_z_N = -(B_inverse_matrix * N_matrix)_transpose * e_i
        std::vector<double> delta_z_N;
        for(int i = 0; i < model.get_num_cols(); ++i){
            double value = 0; 
            for(int j = 0; j < model.get_num_rows(); ++j)
                value += model.get_basis()(leaving_index, j) * model.get_N()(j, i);
            delta_z_N.push_back(-value);
        }
        std::cout<<"finish step 3"<<std::endl;

        // step 4. compute dual step length
        // s = (max(delta_z_j / z_star_j))^-1
        double s = INT_MIN;
        int entering_index = -1;
        for(int i = 0; i < model.get_num_cols(); ++i){
            double ratio = delta_z_N[i]/model.get_header()(i);
            if (ratio > s && ratio != 0){
                s = ratio;
                entering_index = i;
            }
        }
        std::cout<<"finish step 4"<<std::endl;

        // step 5. select entering variable
        if (entering_index == -1) return false;
        s = 1/s;
        SOPVariables entering_variable = model.get_nonbasic(entering_index);
        std::cout<<"finish step 5"<<std::endl;

        // step 6. compute primal step direction
        // delta_x_B = B_matrix * N_matrix * e_j
        std::vector<double> delta_x_B;
        for(int i = 0 ; i < model.get_num_rows(); ++i){
            double value = 0;
            for (int j = 0; j < model.get_num_rows(); ++j)
                value += model.get_basis()(i,j) * model.get_N()(j, entering_index);
            delta_x_B.push_back(value);
        }
        std::cout<<"finish step 6"<<std::endl;

        // step 7. compute primal step length
        double t = model.get_rhs()(leaving_index)/delta_x_B[leaving_index];
        std::cout<<"finish step 7"<<std::endl;

        // step 8. update current primal and dual solution
        for (int i = 0; i < model.get_num_rows(); ++i){
            model.set_rhs(i, model.get_rhs()(i)-t*delta_x_B[i]);
            model.set_header(i, model.get_header()(i)-s*delta_z_N[i]);
        }
        std::cout<<"finish step8"<<std::endl;

        // step 9. update basis
        // use eta_matrix to update B_matrix
        // eta_inverse_matrix = I -((delta_x_B - e_i)e_i_T) / delta_x_i
        // new_B_inverse_matrix = eta_inverse_matrix * old_B_inverse_matrix
        model.set_eta(leaving_index, delta_x_B);
        std::cout<<"finish set_eta function"<<std::endl;
        model.update_basis();
        std::cout<<"finish update basis"<<std::endl;

        // update N_matrix
        for (int i = 0; i < model.get_num_rows(); ++i) 
            model.get_N()(i,entering_index) = model.get_all_range()(i,leaving_variable.get_var_id());
        std::cout<<"finish N matrix"<<std::endl;

        // update nonbasic and basic set
        model.set_basic(leaving_index, entering_variable);
        model.set_nonbasic(entering_index, leaving_variable);

        model.set_rhs(leaving_index, t);
        model.set_header(entering_index, s);
        std::cout<< "entering variable index:" << entering_index << ", name: " << entering_variable.get_var_name() << std::endl;
        std::cout<< "leaving variable index:" << leaving_index << ", name: " << leaving_variable.get_var_name() << std::endl;
        std::cout<<"finish step 9"<<std::endl;

        this->d_iterations++;
    }
    return true;
}

bool SOPlpSolve::netword_simplex(){
    return true;
}

bool SOPlpSolve::barrier_method(){
    return true;
}

bool SOPlpSolve::solve(SOPModel &model, int algorithm){
    bool res;
    switch (algorithm)
    {
    case 1:
        res = primal_simplex(model);
        break;
    case 2:
        res = dual_simplex(model);
        break;
    case 3:
        res = netword_simplex();
        break;
    case 4:
        res = barrier_method();
        break;
    default:
        res = primal_simplex(model);
        break;
    }
    return res;
}
