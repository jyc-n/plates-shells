#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
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

void config_bc(const Parameters& SimPar, std::vector<bool>& bc_d, Eigen::VectorXd& bc_f);

void calculate_fstretch(const Parameters& SimPar, Element& el, Eigen::VectorXd& dEdq);

void calculate_fshear(const Parameters& SimPar, Element& el, Eigen::VectorXd& dEdq);

void calculate_fbend(const Parameters& SimPar, Element& el, Eigen::VectorXd& dEdq);

void calculate_jstretch(const Parameters& SimPar, Element& el, Eigen::MatrixXd& ddEddq);

void calculate_jshear(const Parameters& SimPar, Element& el, Eigen::MatrixXd& ddEddq);

void calculate_jbend(const Parameters& SimPar, Element& el, Eigen::MatrixXd& ddEddq);

void update_residual(const Parameters& SimPar,
                   Eigen::VectorXd& qn,
                   Eigen::VectorXd& qnew,
                   Eigen::VectorXd& vel,
                   Eigen::VectorXd& ext_f,
                   Eigen::VectorXd& res_f,
                   Element* l_element);

void update_jacobian(const Parameters& SimPar,
                   Eigen::VectorXd& qn,
                   Eigen::VectorXd& qnew,
                   Eigen::MatrixXd& mat_j,
                   Element* l_element);

void calculate_dofnew(const Parameters& SimPar, Eigen::VectorXd& qn, Eigen::VectorXd& qnew,
                      std::vector<bool>& bc_d, Eigen::VectorXd& res_f, Eigen::MatrixXd& mat_j);

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
    Eigen::VectorXd vel = Eigen::VectorXd::Zero(SimPar.nn() * SimPar.ndof());

    // boundary conditions
    Eigen::VectorXd fext = Eigen::VectorXd::Zero(SimPar.nn() * SimPar.ndof());
    std::vector<bool> disp(SimPar.nn() * SimPar.ndof(), true);
    config_bc(SimPar, disp, fext);

    // initial guess q_{n+1} = q_n
    dof = m_dof;
    dof_new = dof;


    // TODO: implement a timer class

    auto t_start = system_clock::now();
    // calculate initial residual vector
    update_residual(SimPar, dof, dof_new, vel, fext, residual, m_element_lst);
    auto t_end = system_clock::now();
    auto t_duration = duration_cast<milliseconds>(t_end - t_start);
    std::cout << "residual calculation finished in " << double(t_duration.count()) << " ms" << std::endl;
    std::cout << "residual is " << residual.norm() << std::endl;

    // calculate initial jacobian matrix
    update_jacobian(SimPar, dof, dof_new, jacobian, m_element_lst);

    // tolerance
    double tol = SimPar.kstretch() * 1e-4;

    for (int i = 0; i < SimPar.nst(); i++) {
        std::cout << "--------Step " << i << "--------" << std::endl;

        // convergence flag
        bool CONVERGED = false;

        // apply Newton-Raphson Method
        for (int niter = 0; niter < SimPar.iter_lim(); niter++) {

            // convergence criteria
            if (residual.norm() < tol) {
                std::cout << "Newton's method converges in " << niter+1 << " step(s)" << std::endl;
                CONVERGED = true;
                break;
            }

            t_start = system_clock::now();
            // calculate residual vector
            update_residual(SimPar, dof, dof_new, vel, fext, residual, m_element_lst);
            t_end = system_clock::now();
            t_duration = duration_cast<milliseconds>(t_end - t_start);
            std::cout << "residual calculation finished in " << double(t_duration.count()) << " ms" << std::endl;
            std::cout << "residual is " << residual.norm() << std::endl;

            // calculate jacobian matrix
            update_jacobian(SimPar, dof, dof_new, jacobian, m_element_lst);

            // solve for new dof vector
            calculate_dofnew(SimPar, dof, dof_new, disp, residual, jacobian);
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

// set up boundary conditions
void config_bc(const Parameters& SimPar, std::vector<bool>& bc_d, Eigen::VectorXd& bc_f) {

    // fixed left edge
    for (unsigned int nn = 1; nn <= SimPar.nn() - SimPar.num_nodes_len() + 1; nn += SimPar.num_nodes_len()) {
        bc_d[3*(nn-1)] = false;
        bc_d[3*(nn-1)+1] = false;
        bc_d[3*(nn-1)+2] = false;
    }

    // apply force on the right edge, Fy = 10
    for (unsigned int nn = SimPar.num_nodes_len(); nn <= SimPar.nn(); nn += SimPar.num_nodes_len()) {
        bc_f(3*(nn-1)+1) = 10;
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
                break;
            case 2:
                loc1 = 2;
                loc2 = 3;
                break;
            case 3:
                loc1 = 1;
                loc2 = 3;
                break;
        }

        // coordinates of local nodes
        Eigen::Vector3d loc_xyz1 = el.get_node(loc1)->get_xyz();
        Eigen::Vector3d loc_xyz2 = el.get_node(loc2)->get_xyz();

        // local node number corresponds to global node number
        unsigned int n1 = el.get_node_num(loc1);
        unsigned int n2 = el.get_node_num(loc2);

        // stretch force for local node1
        dEdq(3*(n1-1))   += dEs_dx1(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        dEdq(3*(n1-1)+1) += dEs_dy1(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        dEdq(3*(n1-1)+2) += dEs_dz1(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());

        // stretch force for local node2
        dEdq(3*(n2-1))   += dEs_dx2(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        dEdq(3*(n2-1)+1) += dEs_dy2(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        dEdq(3*(n2-1)+2) += dEs_dz2(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
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
    dEdq(3*(n1-1))   += dEsh_dx1(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kshear());
    dEdq(3*(n1-1)+1) += dEsh_dy1(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kshear());
    dEdq(3*(n1-1)+2) += dEsh_dz1(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kshear());

    // stretch force for local node2
    dEdq(3*(n2-1))   += dEsh_dx2(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kshear());
    dEdq(3*(n2-1)+1) += dEsh_dy2(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kshear());
    dEdq(3*(n2-1)+2) += dEsh_dz2(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kshear());

    // stretch force for local node3
    dEdq(3*(n3-1))   += dEsh_dx3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kshear());
    dEdq(3*(n3-1)+1) += dEsh_dy3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kshear());
    dEdq(3*(n3-1)+2) += dEsh_dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kshear());

}

// TODO: modify bending derivatives and implement in this code
void calculate_fbend(const Parameters& SimPar, Element& el, Eigen::VectorXd& dEdq);

void calculate_jstretch(const Parameters& SimPar, Element& el, Eigen::MatrixXd& ddEddq) {

    // loop over 3 edges
    for (int j = 1; j <= 3; j++) {

        double l0 = el.get_len_edge(j);
        int loc1, loc2 = 0;

        switch (j) {
            case 1:
                loc1 = 1;
                loc2 = 2;
                break;
            case 2:
                loc1 = 2;
                loc2 = 3;
                break;
            case 3:
                loc1 = 1;
                loc2 = 3;
                break;
        }

        // coordinates of local nodes
        Eigen::Vector3d loc_xyz1 = el.get_node(loc1)->get_xyz();
        Eigen::Vector3d loc_xyz2 = el.get_node(loc2)->get_xyz();

        // local jacobian matrix
        Eigen::MatrixXd loc_j = Eigen::MatrixXd::Zero(2 * SimPar.ndof(), 2 * SimPar.ndof());

        loc_j(0, 0) = ddEs_dx1dx1(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        loc_j(0, 1) = ddEs_dx1dy1(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        loc_j(0, 2) = ddEs_dx1dz1(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        loc_j(0, 3) = ddEs_dx1dx2(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        loc_j(0, 4) = ddEs_dx1dy2(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        loc_j(0, 5) = ddEs_dx1dz2(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());

        loc_j(1, 1) = ddEs_dy1dy1(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        loc_j(1, 2) = ddEs_dy1dz1(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        loc_j(1, 3) = ddEs_dx2dy1(loc_xyz1, loc_xyz2, l0, SimPar.kstretch()); //
        loc_j(1, 4) = ddEs_dy1dy2(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        loc_j(1, 5) = ddEs_dy1dz2(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());

        loc_j(2, 2) = ddEs_dz1dz1(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        loc_j(2, 3) = ddEs_dx2dz1(loc_xyz1, loc_xyz2, l0, SimPar.kstretch()); //
        loc_j(2, 4) = ddEs_dy2dz1(loc_xyz1, loc_xyz2, l0, SimPar.kstretch()); //
        loc_j(2, 5) = ddEs_dz1dz2(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());

        loc_j(3, 3) = ddEs_dx2dx2(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        loc_j(3, 4) = ddEs_dx2dy2(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        loc_j(3, 5) = ddEs_dx2dz2(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());

        loc_j(4, 4) = ddEs_dy2dy2(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        loc_j(4, 5) = ddEs_dy2dz2(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());
        loc_j(5, 5) = ddEs_dz2dz2(loc_xyz1, loc_xyz2, l0, SimPar.kstretch());

        for (int i = 0; i < 2 * SimPar.ndof(); i++) {
            for (int j = 0; j < 2 * SimPar.ndof(); j++) {
                if (i > j)
                    loc_j(i,j) = loc_j(j,i);
            }
        }


        // local node number corresponds to global node number
        unsigned int n1 = el.get_node_num(loc1);
        unsigned int n2 = el.get_node_num(loc2);

        unsigned int nx1 = 3*(n1-1);
        unsigned int nx2 = 3*(n2-1);

        // update jacobian
        for (int i = 0; i < SimPar.ndof(); i++) {
            for (int j = 0; j < SimPar.ndof(); j++) {
                ddEddq(nx1+i, nx1+j) += loc_j(i, j);
                ddEddq(nx1+i, nx2+j) += loc_j(i, j+SimPar.ndof());
                ddEddq(nx2+i, nx1+j) += loc_j(i+SimPar.ndof(), j);
                ddEddq(nx2+i, nx2+j) += loc_j(i+SimPar.ndof(), j+SimPar.ndof());
            }
        }
    }
}

void calculate_jshear(const Parameters& SimPar, Element& el, Eigen::MatrixXd& ddEddq) {

    // coordinates of local nodes
    Eigen::Vector3d loc_xyz1 = el.get_node(1)->get_xyz();
    Eigen::Vector3d loc_xyz2 = el.get_node(2)->get_xyz();
    Eigen::Vector3d loc_xyz3 = el.get_node(3)->get_xyz();

    // area
    double a0 = el.get_area();

    // local jacobian matrix
    Eigen::MatrixXd loc_j = Eigen::MatrixXd::Zero(3 * SimPar.ndof(), 3 * SimPar.ndof());

    loc_j(0, 0) = ddEsh_dx1dx1(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(0, 1) = ddEsh_dx1dy1(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(0, 2) = ddEsh_dx1dz1(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(0, 3) = ddEsh_dx1dx2(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(0, 4) = ddEsh_dx1dy2(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(0, 5) = ddEsh_dx1dz2(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(0, 6) = ddEsh_dx1dx3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(0, 7) = ddEsh_dx1dy3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(0, 8) = ddEsh_dx1dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());

    loc_j(1, 1) = ddEsh_dy1dy1(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(1, 2) = ddEsh_dy1dz1(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(1, 3) = ddEsh_dx2dy1(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch()); //
    loc_j(1, 4) = ddEsh_dy1dy2(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(1, 5) = ddEsh_dy1dz2(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(1, 6) = ddEsh_dx3dy1(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch()); //
    loc_j(1, 7) = ddEsh_dy1dy3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(1, 8) = ddEsh_dy1dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());

    loc_j(2, 2) = ddEsh_dz1dz1(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(2, 3) = ddEsh_dx2dz1(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch()); //
    loc_j(2, 4) = ddEsh_dy2dz1(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch()); //
    loc_j(2, 5) = ddEsh_dz1dz2(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(2, 6) = ddEsh_dx3dz1(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch()); //
    loc_j(2, 7) = ddEsh_dy3dz1(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch()); //
    loc_j(2, 8) = ddEsh_dz1dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());

    loc_j(3, 3) = ddEsh_dx2dx2(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(3, 4) = ddEsh_dx2dy2(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(3, 5) = ddEsh_dx2dz2(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(3, 6) = ddEsh_dx2dx3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(3, 7) = ddEsh_dx2dy3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(3, 8) = ddEsh_dx2dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());

    loc_j(4, 4) = ddEsh_dy2dy2(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(4, 5) = ddEsh_dy2dz2(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(4, 6) = ddEsh_dx2dx3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(4, 7) = ddEsh_dx2dy3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(4, 8) = ddEsh_dx2dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());

    loc_j(5, 5) = ddEsh_dz2dz2(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(5, 6) = ddEsh_dx3dz2(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch()); //
    loc_j(5, 7) = ddEsh_dy3dz2(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch()); //
    loc_j(5, 8) = ddEsh_dz2dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());

    loc_j(6, 6) = ddEsh_dx3dx3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(6, 7) = ddEsh_dx3dy3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(6, 8) = ddEsh_dx3dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());

    loc_j(7, 7) = ddEsh_dy3dy3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(7, 8) = ddEsh_dy3dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());
    loc_j(8, 8) = ddEsh_dz3dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, SimPar.kstretch());

    for (int i = 0; i < 3 * SimPar.ndof(); i++) {
        for (int j = 0; j < 3 * SimPar.ndof(); j++) {
            if (i > j)
                loc_j(i,j) = loc_j(j,i);
        }
    }

    // local node number corresponds to global node number
    unsigned int n1 = el.get_node_num(1);
    unsigned int n2 = el.get_node_num(2);
    unsigned int n3 = el.get_node_num(3);

    unsigned int nx1 = 3*(n1-1);
    unsigned int nx2 = 3*(n2-1);
    unsigned int nx3 = 3*(n3-1);

    // update jacobian
    for (int i = 0; i < SimPar.ndof(); i++) {
        for (int j = 0; j < SimPar.ndof(); j++) {

            ddEddq(nx1+i, nx1+j) += loc_j(i, j);
            ddEddq(nx1+i, nx2+j) += loc_j(i, j+SimPar.ndof());
            ddEddq(nx1+i, nx3+j) += loc_j(i, j+SimPar.ndof()*2);

            ddEddq(nx2+i, nx1+j) += loc_j(i+SimPar.ndof(), j);
            ddEddq(nx2+i, nx2+j) += loc_j(i+SimPar.ndof(), j+SimPar.ndof());
            ddEddq(nx2+i, nx3+j) += loc_j(i+SimPar.ndof(), j+SimPar.ndof()*2);

            ddEddq(nx3+i, nx1+j) += loc_j(i+SimPar.ndof()*2, j);
            ddEddq(nx3+i, nx2+j) += loc_j(i+SimPar.ndof()*2, j+SimPar.ndof());
            ddEddq(nx3+i, nx3+j) += loc_j(i+SimPar.ndof()*2, j+SimPar.ndof()*2);
        }
    }
}

void calculate_jbend(const Parameters& SimPar, Element& el, Eigen::VectorXd& ddEddq);

/*
 *    f_i = m_i * (q_i(t_n+1) - q_i(t_n)) / dt^2 - m_i * v(t_n) / dt + dE/dt - F_ext
 */
void update_residual(const Parameters& SimPar,
                     Eigen::VectorXd& qn,
                     Eigen::VectorXd& qnew,
                     Eigen::VectorXd& vel,
                     Eigen::VectorXd& ext_f,
                     Eigen::VectorXd& res_f,
                     Element* l_element) {

    // 1st derivative of energy vector
    Eigen::VectorXd dEdq = Eigen::VectorXd::Zero(SimPar.nn() * SimPar.ndof());

    // loop over each element
    for (unsigned int i = 0; i < SimPar.nel(); i++) {

        calculate_fstretch(SimPar, l_element[i], dEdq);
        calculate_fshear(SimPar, l_element[i], dEdq);
        // TODO: finish bending
    }

    double dt = SimPar.dt();

    // TODO: calculate mass
    double mi = 1.0 / 6.0;

    for (unsigned int i = 0; i < SimPar.nn() * SimPar.ndof(); i++)
    {
        res_f(i) = mi * (qnew(i) - qn(i)) / pow(dt,2) - mi * vel(i)/dt + dEdq(i) - ext_f(i);
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

    // 2nd derivative of energy vector
    Eigen::MatrixXd ddEddq = Eigen::MatrixXd::Zero(SimPar.nn() * SimPar.ndof(), SimPar.nn() * SimPar.ndof());
    for (unsigned int i = 0; i < SimPar.nel(); i++) {

        calculate_jstretch(SimPar, l_element[i], ddEddq);
        //calculate_jshear(SimPar, l_element[i], ddEddq);
        // TODO: finish bending
    }
}

/*
 *    q_{n+1} = q_n - J \ f
 */
void calculate_dofnew(const Parameters& SimPar, Eigen::VectorXd& qn, Eigen::VectorXd& qnew,
                      std::vector<bool>& bc_d, Eigen::VectorXd& res_f, Eigen::MatrixXd& mat_j) {

    int num_fixed = 0;
    int num_total = SimPar.nn() * SimPar.ndof();
    int num_free = num_total - num_fixed;

    std::vector<int> fix_dof;
    for (int i = 0; i < bc_d.size(); i++) {
         if ( !bc_d[i] ) {
             num_fixed++;
             fix_dof.push_back(i);
         }
    }

    Eigen::VectorXd dq = Eigen::VectorXd::Zero(num_free);
    Eigen::VectorXd temp_f = Eigen::VectorXd::Zero(num_free);
    Eigen::MatrixXd temp_j = Eigen::MatrixXd::Zero(num_free, num_free);

    for (int i = 0, j = 0, k = 0; i < num_free && j < num_total && k < fix_dof.size(); j++) {
        if (fix_dof[k] == j) {
            k++;
            continue;
        }
        temp_f(i) = res_f(j);
    }

    for (int i1 = 0, j1 = 0, k = 0; i1 < num_free  && j1 < num_total && k < fix_dof.size(); j1++) {
        for (int i2 = 0, j2 = 0; i2 < num_free && j2 < num_total; j2++) {

            if (fix_dof[k] == j1 || fix_dof[k] == j2) {
                k++;
                continue;
            }
            temp_f(i1, i2) = mat_j(j1, j2);
        }
    }

    dq = temp_j.householderQr().solve(temp_f);

    std::cout << dq.norm() << std::endl;


}