#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include "solver.h"
#include "element.h"
#include "geometry.h"

// ========================================= //
//         Declaration of subroutines        //
// ========================================= //

void prepare_solver(const Parameters& SimPar,
                    Node* l_node,
                    Element* l_element);

void update_residual(const Parameters& SimPar,
                   Eigen::VectorXd& qn,
                   Eigen::VectorXd& qnew,
                   Eigen::VectorXd& vel,
                   Eigen::VectorXd& res_f);

void update_jacobian(const Parameters& SimPar,
                   Eigen::VectorXd& qn,
                   Eigen::VectorXd& qnew,
                   Eigen::MatrixXd& mat_j);

// ========================================= //
//            Main solver part               //
// ========================================= //

void fake_solver(Parameters& SimPar) {
    /*
    double ti = 0.0;
    double tf = 2.0;
    double dt = 0.1;
    double amp = 2.0 * dt / (tf - ti);

    Eigen::MatrixXd fake_motion = Eigen::MatrixXd::Zero(SimPar.nn(), SimPar.ndof());
    for (int i = 0; i < SimPar.nn(); i++) {
        for (int j = 0; j < SimPar.ndof(); j++) {
            if (j == 2) { // z-direction
                fake_motion(i,j) = amp * sin(M_PI / 4.0 * InitGeo.lst_coord()(i,0));
            }
        }
    }
    std::cout << fake_motion << '\n';

    Eigen::MatrixXd coord = InitGeo.lst_coord();
    int counter = 0;
    for (double t = ti; t < tf + 0.0000001; t += dt) {
        std::string filepath = "/Users/chenjingyu/git/plates-shells/";
        std::string filename = "result" + std::to_string(counter) + ".txt";
        std::ofstream myfile((filepath+filename).c_str());
        for (int i = 0; i < SimPar.nn(); i++) {
            for (int j = 0; j < SimPar.ndof(); j++) {
                myfile << std::setprecision(6) << std::fixed << coord(i,j) << '\t';
            }
            myfile << std::endl;
        }
        coord += fake_motion;
        counter++;
        //std::cout << '\n';
    }*/
}

void time_solver(Parameters& SimPar) {
    Node m_node_lst[SimPar.nn()];
    Element m_element_lst[SimPar.nel()];

    prepare_solver(SimPar, m_node_lst, m_element_lst);

    // solver starts here

    // residual vector
    Eigen::VectorXd residual = Eigen::VectorXd::Zero(SimPar.nn() * SimPar.ndof());

    // Jacobian matrix
    Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(SimPar.nn() * SimPar.ndof(), SimPar.nn() * SimPar.ndof());

    // dof vector at t_n and t_{n+1}
    Eigen::VectorXd dof(SimPar.nn() * SimPar.ndof());
    Eigen::VectorXd dof_new(SimPar.nn() * SimPar.ndof());

    // velocity vector
    Eigen::VectorXd vel(SimPar.nn() * SimPar.ndof());

    // initial guess q_{n+1} = q_n
    dof = m_dof;
    dof_new = dof;

    // calculate initial residual vector
    update_residual(SimPar, dof, dof_new, vel, residual);

    // calculate initial jacobian matrix
    update_jacobian(SimPar, dof, dof_new, jacobian);

    for (int i = 0; i < SimPar.nst(); i++) {
        std::cout << "--------Step " << i << "--------" << std::endl;

        // convergence flag
        bool CONVERGED = false;

        // apply Newton-Raphson Method
        for (int niter = 0; niter < SimPar.iter_lim(); niter++) {

            // convergence criteria
            // TODO: add a parameter for convergence criteria
            if (residual.norm() < 1e-6) {
                std::cout << "Newton's method converges in " << niter+1 << " step(s)" << std::endl;
                CONVERGED = true;
                break;
            }

            // calculate residual vector
            update_residual(SimPar, dof, dof_new, vel, residual);

            // calculate jacobian matrix
            update_jacobian(SimPar, dof, dof_new, jacobian);

            /*
             *    q_{n+1} = q_n - J \ f
             */

            std::cout << "pretend doing sth" << std::endl;
        }
        if (!CONVERGED) {
            std::cerr << "Solver did not converge in " << SimPar.iter_lim()
                      << " iterations at step " << i << std::endl;
            throw "Cannot converge! Program terminated";
        }
    }

}

// ========================================= //
//       Implementation of subroutines       //
// ========================================= //

// initialize element list
void prepare_solver(const Parameters& SimPar, Node* l_node, Element* l_element) {
    for (unsigned int i = 0; i < SimPar.nn(); i++) {
        Node temp(i+1, m_coord(i,0), m_coord(i,1), m_coord(i,2));
        l_node[i] = temp;
    }

    for (unsigned int i = 0; i < SimPar.nel(); i++) {
        // pointers to 3 node object
        Node* pn1 = &l_node[m_conn(i,0)-1];
        Node* pn2 = &l_node[m_conn(i,1)-1];
        Node* pn3 = &l_node[m_conn(i,2)-1];

        Element temp(i+1, pn1, pn2, pn3);

        temp.find_nearby_element(SimPar);
        l_element[i] = temp;
    }
}


/*
 *    f_i = m_i * (q_i(t_n+1) - q_i(t_n)) / dt^2 - m_i * v(t_n) / dt + dE/dt - F_ext
 */
void update_residual(const Parameters& SimPar,
                     Eigen::VectorXd& qn,
                     Eigen::VectorXd& qnew,
                     Eigen::VectorXd& vel,
                     Eigen::VectorXd& res_f) {
    //TODO: implement the equation 16 and global assembly
}



/*
 *    J_ij = m_i / dt^2 * delta_ij + d^2 E / dq_i dq_j
 */
void update_jacobian(const Parameters& SimPar,
                     Eigen::VectorXd& qn,
                     Eigen::VectorXd& qnew,
                     Eigen::MatrixXd& mat_j) {
    //TODO: implement the equation 17 and global assembly
}