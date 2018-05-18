#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <chrono>
#include <Eigen/Dense>
#include "solver.h"
#include "element.h"
#include "geometry.h"
#include "stretch_derivatives.h"
#include "shear_derivatives.h"

using namespace std::chrono;

// ========================================= //
//         Declaration of subroutines        //
// ========================================= //

void prepare_solver(const Parameters& SimPar, Node* l_node, Element* l_element);

void calculate_fstretch(const Parameters& SimPar, Element& el, Eigen::VectorXd& dEdq);

void calculate_fshear(const Parameters& SimPar, Element& el, Eigen::VectorXd& dEdq);

void calculate_fbend(const Parameters& SimPar, Element& el, Eigen::VectorXd& dEdq);

void calculate_jstretch(const Parameters& SimPar, Element& el, Eigen::VectorXd& ddEddq);

void calculate_jshear(const Parameters& SimPar, Element& el, Eigen::VectorXd& ddEddq);

void calculate_jbend(const Parameters& SimPar, Element& el, Eigen::VectorXd& ddEddq);

void update_residual(const Parameters& SimPar,
                   Eigen::VectorXd& qn,
                   Eigen::VectorXd& qnew,
                   Eigen::VectorXd& vel,
                   Eigen::VectorXd& res_f,
                   Element* l_element);

void update_jacobian(const Parameters& SimPar,
                   Eigen::VectorXd& qn,
                   Eigen::VectorXd& qnew,
                   Eigen::MatrixXd& mat_j,
                   Element* l_element);

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

    auto t_start = system_clock::now();
    // calculate initial residual vector
    update_residual(SimPar, dof, dof_new, vel, residual, m_element_lst);
    auto t_end = system_clock::now();
    auto t_duration = duration_cast<microseconds>(t_end - t_start);
    std::cout << "residual calculation finished in " << double(t_duration.count()) << " ms" << std::endl;

    // calculate initial jacobian matrix
    update_jacobian(SimPar, dof, dof_new, jacobian, m_element_lst);

    for (int i = 0; i < SimPar.nst(); i++) {
        std::cout << "--------Step " << i << "--------" << std::endl;

        // convergence flag
        bool CONVERGED = false;

        // apply Newton-Raphson Method
        for (int niter = 0; niter < SimPar.iter_lim(); niter++) {

            // convergence criteria
            // TODO: implement dynamic convergence criteria

            if (residual.norm() < 1e-6) {
                std::cout << "Newton's method converges in " << niter+1 << " step(s)" << std::endl;
                CONVERGED = true;
                break;
            }

            auto t_start = system_clock::now();
            // calculate residual vector
            update_residual(SimPar, dof, dof_new, vel, residual, m_element_lst);
            auto t_end = system_clock::now();
            auto t_duration = duration_cast<microseconds>(t_end - t_start);
            std::cout << "residual calculation finished in " << double(t_duration.count()) << " ms" << std::endl;

            // calculate jacobian matrix
            update_jacobian(SimPar, dof, dof_new, jacobian, m_element_lst);

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

// stretch energy for each element
void calculate_fstretch(const Parameters& SimPar, Element& el, Eigen::VectorXd& dEdq) {

    // loop over 3 edges
    for (int j = 1; j <= 3; j++) {

        double l0 = el.get_len_edge(j);
        int loc1, loc2 = 0;

        switch (j) {
            case 1:
                loc1 = 1;
                loc2 = 2;
            case 2:
                loc1 = 2;
                loc2 = 3;
            case 3:
                loc1 = 1;
                loc2 = 3;
        }

        // coordinates of local nodes
        Eigen::Vector3d loc_xyz1 = el.get_node(loc1)->get_xyz();
        Eigen::Vector3d loc_xyz2 = el.get_node(loc2)->get_xyz();

        // local node number corresponds to global node number
        unsigned int n1 = el.get_node_num(loc1);
        unsigned int n2 = el.get_node_num(loc2);

        // stretch force for local node1
        dEdq[3*(n1-1)]   += dEs_dx1(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        dEdq[3*(n1-1)+1] += dEs_dy1(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        dEdq[3*(n1-1)+2] += dEs_dz1(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());

        // stretch force for local node2
        dEdq[3*(n2-1)]   += dEs_dx2(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        dEdq[3*(n2-1)+1] += dEs_dy2(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        dEdq[3*(n2-1)+2] += dEs_dz2(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
    }
}

// shear energy for each element
void calculate_fshear(const Parameters& SimPar, Element& el, Eigen::VectorXd& dEdq) {

    // coordinates of local nodes
    Eigen::Vector3d loc_xyz1 = el.get_node(1)->get_xyz();
    Eigen::Vector3d loc_xyz2 = el.get_node(2)->get_xyz();
    Eigen::Vector3d loc_xyz3 = el.get_node(3)->get_xyz();

    // local node number corresponds to global node number
    unsigned int n1 = el.get_node_num(1);
    unsigned int n2 = el.get_node_num(2);
    unsigned int n3 = el.get_node_num(3);

    // area
    double a0 = el.get_area();

    // stretch force for local node1
    dEdq[3*(n1-1)]   += dEsh_dx1(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kshear());
    dEdq[3*(n1-1)+1] += dEsh_dy1(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kshear());
    dEdq[3*(n1-1)+2] += dEsh_dz1(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kshear());

    // stretch force for local node2
    dEdq[3*(n2-1)]   += dEsh_dx2(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kshear());
    dEdq[3*(n2-1)+1] += dEsh_dy2(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kshear());
    dEdq[3*(n2-1)+2] += dEsh_dz2(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kshear());

    // stretch force for local node3
    dEdq[3*(n3-1)]   += dEsh_dx3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kshear());
    dEdq[3*(n3-1)+1] += dEsh_dy3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kshear());
    dEdq[3*(n3-1)+2] += dEsh_dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kshear());

}

void calculate_fbend(const Parameters& SimPar, Element& el, Eigen::VectorXd& dEdq);

void calculate_jstretch(const Parameters& SimPar, Element& el, Eigen::VectorXd& ddEddq);

void calculate_jshear(const Parameters& SimPar, Element& el, Eigen::VectorXd& ddEddq);

void calculate_jbend(const Parameters& SimPar, Element& el, Eigen::VectorXd& ddEddq);

/*
 *    f_i = m_i * (q_i(t_n+1) - q_i(t_n)) / dt^2 - m_i * v(t_n) / dt + dE/dt - F_ext
 */
void update_residual(const Parameters& SimPar,
                     Eigen::VectorXd& qn,
                     Eigen::VectorXd& qnew,
                     Eigen::VectorXd& vel,
                     Eigen::VectorXd& res_f,
                     Element* l_element) {
    //TODO: implement the equation 16 and global assembly

    // 1st derivative of energy vector
    Eigen::VectorXd dEdq = Eigen::VectorXd::Zero(SimPar.nn() * SimPar.ndof());

    // loop over each element
    for (unsigned int i = 0; i < SimPar.nel(); i++) {

        // TODO: initial length/area/angle should be passed in
        calculate_fstretch(SimPar, l_element[i], dEdq);
        calculate_fshear(SimPar, l_element[i], dEdq);
    }

}



/*
 *    J_ij = m_i / dt^2 * delta_ij + d^2 E / dq_i dq_j
 */
void update_jacobian(const Parameters& SimPar,
                     Eigen::VectorXd& qn,
                     Eigen::VectorXd& qnew,
                     Eigen::MatrixXd& mat_j,
                     Element* l_element) {
    //TODO: implement the equation 17 and global assembly
}