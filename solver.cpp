#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <Eigen/Dense>

#include "solver.h"
#include "parameters.h"
#include "geometry.h"
#include "loadbc.h"
#include "node.h"
#include "hinge.h"
#include "element.h"
#include "utilities.h"
#include "stretch_derivatives.h"
#include "shear_derivatives.h"
#include "bending_derivatives.h"


SolverImpl::SolverImpl(Parameters* SimPar, Geometry* SimGeo, Boundary* SimBC) {
    m_SimPar = SimPar;
    m_SimGeo = SimGeo;
    m_SimBC  = SimBC;
}

// ========================================= //
//           Main solver function            //
// ========================================= //

void SolverImpl::Solve() {

    // dof vector at t_n and t_{n+1}
    Eigen::VectorXd dof(m_SimGeo->nn() * m_SimGeo->ndof());
    Eigen::VectorXd dof_new(m_SimGeo->nn() * m_SimGeo->ndof());
    dof = m_SimGeo->m_dof;

    // additional dof vector for Newton's method
    Eigen::VectorXd dof_old(m_SimGeo->nn() * m_SimGeo->ndof());

    // velocity vector
    Eigen::VectorXd vel = Eigen::VectorXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof());

    Timer t;
    Timer t_all(true);

    // TODO: implement tolerance in parameter
    // tolerance
    //double tol = SimPar.kstretch() * 1e-2;
    double tol = 9.81 * 1e-1;

    // true - write output, false - no output, configured in input.txt
    bool WRITE_OUTPUT = m_SimPar->outop();

    for (int i = 0; i < m_SimPar->nst(); i++) {

        // initial guess
        dof_new = dof;
        dof_old = dof;

        std::cout << "--------Step " << i << "--------" << std::endl;
        std::cout << "#iter " << '\t' << "norm of residual" << '\t' << "time" << std::endl;

        // convergence flag
        bool CONVERGED = false;

        // apply Newton-Raphson Method
        for (int niter = 0; niter < m_SimPar->iter_lim(); niter++) {

            // residual vector
            Eigen::VectorXd residual = Eigen::VectorXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof());

            // Jacobian matrix
            Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof(), m_SimGeo->nn() * m_SimGeo->ndof());

            std::cout << niter+1 << '\t';

            t.start();

            Eigen::VectorXd dof_free = unconsVec(dof);
            Eigen::VectorXd dof_new_free = unconsVec(dof_new);
            //std::cout << "dof_free" << std::endl;
            //std::cout << dof_free << std::endl;
            //std::cout << "dof_new_free" << std::endl;
            //std::cout << dof_new_free << std::endl;

            // calculate residual vector
            updateResidual(dof, dof_new, vel, residual);
            //std::cout << residual << std::endl;

            // calculate jacobian matrix
            updateJacobian(jacobian);
            //std::cout << jacobian << std::endl;

            // residual vector and jacobian matrix at dofs that are not specified
            Eigen::VectorXd res_free = unconsVec(residual);
            Eigen::MatrixXd jacobian_free = unconsMat(jacobian);
            //std::cout << res_free << std::endl;

            // solve for new dof vector
            dof_new = calcDofnew(dof_old, res_free, jacobian_free);

            // update coordinates for the nodes
            updateNodes(dof_new);

            // update for next iteration
            dof_old = dof_new;

            // display residual
            std::cout << res_free.norm() << '\t' << t.elapsed() << " ms" << std::endl;

            // convergence criteria
            if (res_free.norm() < tol) {
                std::cout << "Newton's method converges in " << niter + 1 << " step(s)" << std::endl;
                CONVERGED = true;
                break;
            }
        }

        // calculate new velocity vector
        vel = calcVel(m_SimPar->dt(), dof, dof_new);

        dof = dof_new;

        if (WRITE_OUTPUT) {
            // output files
            std::string filepath = "/Users/chenjingyu/Dropbox/Research/Codes/plates-shells/results/";
            std::string filename = "result" + std::to_string(i) + ".txt";
            std::ofstream myfile((filepath+filename).c_str());
            for (int k = 0; k < m_SimGeo->nn(); k++) {
                myfile << std::setprecision(8) << std::fixed
                       << dof(3*k) << '\t'
                       << dof(3*k+1) << '\t'
                       << dof(3*k+2) << std::endl;
            }
        }

        if (!CONVERGED) {
            std::cerr << "Solver did not converge in " << m_SimPar->iter_lim()
                      << " iterations at step " << i << std::endl;
            throw "Cannot converge! Program terminated";
        }
    }
    std::cout << "---------------------------" << std::endl;
    std::cout << "Simulation completed" << std::endl;
    std::cout << "Total time used " << t_all.elapsed(true) << " seconds" <<  std::endl;
}

// ========================================= //
//       Implementation of subroutines       //
// ========================================= //

/*
 *    f_i = m_i * (q_i(t_n+1) - q_i(t_n)) / dt^2 - m_i * v(t_n) / dt + dE/dq - F_ext
 */
void SolverImpl::updateResidual(Eigen::VectorXd& qn, Eigen::VectorXd& qnew, Eigen::VectorXd& vel, Eigen::VectorXd& res_f) {

    // 1st derivative of energy vector
    Eigen::VectorXd dEdq = Eigen::VectorXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof());

    Eigen::VectorXd fstretch = Eigen::VectorXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof());
    Eigen::VectorXd fshear   = Eigen::VectorXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof());
    Eigen::VectorXd fbend    = Eigen::VectorXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof());

    // loop over each element
    for (unsigned int i = 0; i < m_SimGeo->nel(); i++) {
        calculate_fstretch(m_SimGeo->m_elementList[i], fstretch);
        calculate_fshear(m_SimGeo->m_elementList[i], fshear);
        calculate_fbend(m_SimGeo->m_elementList[i], fbend);
    }

/*
    std::cout << "stretch" << std::endl;
    std::cout << fstretch << std::endl;

    std::cout << "shear" << std::endl;
    std::cout << fshear << std::endl;
    std::cout << "bend" << std::endl;
    std::cout << fbend << std::endl;
     */

    dEdq = fstretch + fshear + fbend;


    // reset all marks
    resetMarks();

    double dt = m_SimPar->dt();

    for (unsigned int i = 0; i < m_SimGeo->nn() * m_SimGeo->ndof(); i++)
        res_f(i) = m_SimGeo->m_mass(i) * (qnew(i) - qn(i)) / pow(dt,2) - m_SimGeo->m_mass(i) * vel(i)/dt
                 + dEdq(i) - m_SimBC->m_fext(i);

    /*
    std::cout << "inertial term" << std::endl;
    std::cout << mi * (qnew - qn) / pow(dt,2) << std::endl;
    std::cout << "velocity term" << std::endl;
    std::cout << mi * vel/dt << std::endl;

    std::cout << "dEdq" << std::endl;
    std::cout << dEdq << std::endl;
    std::cout << "residual" << std::endl;
    std::cout << res_f << std::endl;
     */
}



/*
 *    J_ij = m_i / dt^2 * delta_ij + d^2 E / dq_i dq_j
 */
void SolverImpl::updateJacobian(Eigen::MatrixXd& mat_j) {

    // 2nd derivative of energy vector
    Eigen::MatrixXd ddEddq = Eigen::MatrixXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof(), m_SimGeo->nn() * m_SimGeo->ndof());

    Eigen::MatrixXd jstretch = Eigen::MatrixXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof(), m_SimGeo->nn() * m_SimGeo->ndof());
    Eigen::MatrixXd jshear   = Eigen::MatrixXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof(), m_SimGeo->nn() * m_SimGeo->ndof());
    Eigen::MatrixXd jbend    = Eigen::MatrixXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof(), m_SimGeo->nn() * m_SimGeo->ndof());

    for (unsigned int i = 0; i < m_SimGeo->nel(); i++) {
        calculate_jstretch(m_SimGeo->m_elementList[i], jstretch);
        calculate_jshear(m_SimGeo->m_elementList[i], jshear);
        calculate_jbend(m_SimGeo->m_elementList[i], jbend);
    }

    ddEddq = jstretch + jshear + jbend;

/*
    std::cout << "stretch" << std::endl;
    std::cout << jstretch << std::endl;
    std::cout << "shear" << std::endl;
    std::cout << jshear << std::endl;
    std::cout << "bend" << std::endl;
    std::cout << jbend << std::endl;
     */

    // reset all marks
    resetMarks();

    double dt = m_SimPar->dt();

    for (int i = 0; i < m_SimGeo->nn() * m_SimGeo->ndof(); i++) {
        for (int j = 0; j < m_SimGeo->nn() * m_SimGeo->ndof(); j++) {
            mat_j(i, j) = ddEddq(i, j);
        }
        mat_j(i, i) += m_SimGeo->m_mass(i) / pow(dt, 2);
    }
}

Eigen::VectorXd SolverImpl::calcVel(double dt, const Eigen::VectorXd& qcurr, const Eigen::VectorXd& qnew) {
    return (qnew - qcurr) / dt;
}


/*
 *    q_{n+1} = q_n - J \ f
 */
Eigen::VectorXd SolverImpl::calcDofnew(Eigen::VectorXd& qn, const Eigen::VectorXd& temp_f, const Eigen::MatrixXd& temp_j) {

    Eigen::VectorXd dq_free = Eigen::VectorXd::Zero(m_SimBC->m_numFree);
    Eigen::VectorXd dq      = Eigen::VectorXd::Zero(m_SimBC->m_numTotal);

    // 0 * x = 0
    if (abs(temp_f.norm()) < 1e-16 && abs(temp_j.norm()) < 1e-16) {
        std::cout << "Warning: 0 * x = 0 case" << std::endl;
        dq_free.setZero();
    }
        // regular case
    else
        dq_free = temp_j.llt().solve(temp_f);

    //std::cout << res_f << std::endl;

    //std::cout << "---dq_free---" << std::endl;
    //std::cout << dq_free << std::endl;
    //std::cout << "---temp_f---" << std::endl;
    //std::cout << temp_f << std::endl;
    //std::cout << "---temp_j---" << std::endl;
    //std::cout << temp_j << std::endl;

    // map the free part dof vector back to the full dof vector
    for (int i = 0, j = 0; i < m_SimBC->m_numTotal && j < m_SimBC->m_numFree; i++) {
        std::vector<int>::const_iterator it = find(m_SimBC->m_specifiedDof.begin(), m_SimBC->m_specifiedDof.end(), i);
        if (it == m_SimBC->m_specifiedDof.end()) {
            dq(i) = dq_free(j);
            j++;
        }
    }
    //std::cout << dq << std::endl;
    return (qn - dq);
}

// ========================================= //
//          Other helper functions           //
// ========================================= //

Eigen::VectorXd SolverImpl::unconsVec(const Eigen::VectorXd& vec) {

    Eigen::VectorXd vec_free = Eigen::VectorXd::Zero(m_SimBC->m_numFree);

    std::vector<int>::const_iterator p = m_SimBC->m_specifiedDof.begin();
    for (int i = 0, j = 0; i < m_SimBC->m_numTotal; i++) {
        if (i == *p) {
            if (p != m_SimBC->m_specifiedDof.end())
                p++;
            continue;
        }
        else {
            vec_free(j) = vec(i);
            j++;
        }
    }
    return vec_free;
}

Eigen::MatrixXd SolverImpl::unconsMat(const Eigen::MatrixXd& mat) {

    Eigen::MatrixXd mat_free = Eigen::MatrixXd::Zero(m_SimBC->m_numFree, m_SimBC->m_numFree);

    for (int i1 = 0, i2 = 0, j2 = 0; i1 < m_SimBC->m_numTotal; i1++) {
        std::vector<int>::const_iterator pi = find(m_SimBC->m_specifiedDof.begin(), m_SimBC->m_specifiedDof.end(), i1);

        // check if i1 is not fixed
        if (pi == m_SimBC->m_specifiedDof.end()) {

            for (int j1 = 0; j1 < m_SimBC->m_numTotal; j1++) {
                std::vector<int>::const_iterator pj = find(m_SimBC->m_specifiedDof.begin(), m_SimBC->m_specifiedDof.end(), j1);

                // check if j1 is not fixed
                if (pj == m_SimBC->m_specifiedDof.end()) {
                    mat_free(i2, j2) = mat(i1, j1);
                    j2++;
                }
            }
            i2++;
            j2 = 0;
        }
    }
    return mat_free;
}

void SolverImpl::updateNodes(Eigen::VectorXd &qnew){
    for (int i = 0; i < m_SimGeo->nn(); i++) {
        double xpos = qnew(3*i);
        double ypos = qnew(3*i+1);
        double zpos = qnew(3*i+2);
        m_SimGeo->m_nodeList[i].set_xyz(xpos, ypos, zpos);
    }
}

void SolverImpl::resetMarks() {
    for (int i = 0; i < m_SimGeo->nel(); i++) {
        for (int j = 1; j <= 3; j++) {
            // if the j-th edge is a hinge
            if ( m_SimGeo->m_elementList[i].is_hinge(j) )
                m_SimGeo->m_elementList[i].get_hinge(j)->reset_visited();
        }
    }
}

// stretch energy for each element
void SolverImpl::calculate_fstretch(Element& el, Eigen::VectorXd& dEdq) {

    // loop over 3 edges
    for (int i = 1; i <= 3; i++) {

        double l0 = el.get_len_edge(i);
        int loc1 = 0, loc2 = 0;

        // switch over different edges
        switch (i) {
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

        // local stretch force: fx1, fy1, fz1, fx2, fy2, fz2
        Eigen::VectorXd loc_f = Eigen::VectorXd::Zero(6);

        // local node number corresponds to global node number
        unsigned int n1 = el.get_node_num(loc1);
        unsigned int n2 = el.get_node_num(loc2);

        double ks = m_SimPar->kstretch();

        loc_f(0) = dEs_dx1(loc_xyz1, loc_xyz2, l0, ks);
        loc_f(1) = dEs_dy1(loc_xyz1, loc_xyz2, l0, ks);
        loc_f(2) = dEs_dz1(loc_xyz1, loc_xyz2, l0, ks);
        loc_f(3) = dEs_dx2(loc_xyz1, loc_xyz2, l0, ks);
        loc_f(4) = dEs_dy2(loc_xyz1, loc_xyz2, l0, ks);
        loc_f(5) = dEs_dz2(loc_xyz1, loc_xyz2, l0, ks);

        // force needs to be divided by 2, because 1/2 wasn't included during derivation
        loc_f = loc_f * 0.5;

        // if an edge is shared by two elements, local force will be double counted
        if (el.get_nearby_element(i) != 0)
            loc_f = loc_f * 0.5;

        //std::cout << n1 << '\t' << n2 << std::endl;
        //std::cout << loc_f << std::endl;

        // stretch force for local node1
        dEdq(3*(n1-1))   += loc_f(0);
        dEdq(3*(n1-1)+1) += loc_f(1);
        dEdq(3*(n1-1)+2) += loc_f(2);

        // stretch force for local node2
        dEdq(3*(n2-1))   += loc_f(3);
        dEdq(3*(n2-1)+1) += loc_f(4);
        dEdq(3*(n2-1)+2) += loc_f(5);
    }
}

// shear energy for each element
void SolverImpl::calculate_fshear(Element& el, Eigen::VectorXd& dEdq) {

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

    // local shearing force: fx1, fy1, fz1, fx2, fy2, fz2, fx3, fy3, fz3
    Eigen::VectorXd loc_f = Eigen::VectorXd::Zero(9);

    double ksh = m_SimPar->kshear();

    loc_f(0) = dEsh_dx1(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_f(1) = dEsh_dy1(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_f(2) = dEsh_dz1(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_f(3) = dEsh_dx2(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_f(4) = dEsh_dy2(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_f(5) = dEsh_dz2(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_f(6) = dEsh_dx3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_f(7) = dEsh_dy3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_f(8) = dEsh_dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);

    // shearing force for local node1
    dEdq(3*(n1-1))   += loc_f(0);
    dEdq(3*(n1-1)+1) += loc_f(1);
    dEdq(3*(n1-1)+2) += loc_f(2);

    // shearing force for local node2
    dEdq(3*(n2-1))   += loc_f(3);
    dEdq(3*(n2-1)+1) += loc_f(4);
    dEdq(3*(n2-1)+2) += loc_f(5);

    // shearing force for local node3
    dEdq(3*(n3-1))   += loc_f(6);
    dEdq(3*(n3-1)+1) += loc_f(7);
    dEdq(3*(n3-1)+2) += loc_f(8);

}

void SolverImpl::calculate_fbend(Element& el, Eigen::VectorXd& dEdq) {

    // loop over 3 edges
    for (int i = 1; i <= 3; i++) {

        // if the edge is a hinge
        if ( !el.is_hinge(i) )
            continue;

        // get the pointer to that hinge
        Hinge* temp_hinge = el.get_hinge(i);

        // and has NOT been visited yet
        if (temp_hinge->check_visited())
            continue;

        // mark the edge as visited
        temp_hinge->mark_visited();

        // coordinates of local nodes
        Eigen::Vector3d loc_xyz0 = temp_hinge->get_node(0)->get_xyz();
        Eigen::Vector3d loc_xyz1 = temp_hinge->get_node(1)->get_xyz();
        Eigen::Vector3d loc_xyz2 = temp_hinge->get_node(2)->get_xyz();
        Eigen::Vector3d loc_xyz3 = temp_hinge->get_node(3)->get_xyz();

        // local node number corresponds to global node number
        unsigned int n0 = temp_hinge->get_node_num(0);
        unsigned int n1 = temp_hinge->get_node_num(1);
        unsigned int n2 = temp_hinge->get_node_num(2);
        unsigned int n3 = temp_hinge->get_node_num(3);

        // original hinge length
        double l0 = el.get_len_edge(i);

        // original sum of the area of two elements
        double a0 = temp_hinge->get_area_sum();

        // original norm(n1)*norm(n2)
        double n1n20 = temp_hinge->get_n1n20();

        // original kappa0
        Eigen::Vector3d kappa0 = temp_hinge->get_kappa0();

        // local bending force: fx0, fy0, fz0, fx1, fy1, fz1, fx2, fy2, fz2, fx3, fy3, fz3
        Eigen::VectorXd loc_f = Eigen::VectorXd::Zero(12);

        loc_f(0) = dEb_dx0(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_f(1) = dEb_dy0(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_f(2) = dEb_dz0(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_f(3) = dEb_dx1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_f(4) = dEb_dy1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_f(5) = dEb_dz1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_f(6) = dEb_dx2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_f(7) = dEb_dy2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_f(8) = dEb_dz2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_f(9) = dEb_dx3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_f(10) = dEb_dy3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_f(11) = dEb_dz3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);

        // force needs to be multiplied by the coefficients
        loc_f = loc_f * 0.5 * m_SimPar->kbend() * 3 * pow(l0, 2) / a0;

        //std::cout << "bending local force" << std::endl;
        //std::cout << loc_f << std::endl;

        // bending force for local node1
        dEdq(3*(n0-1))   += loc_f(0);
        dEdq(3*(n0-1)+1) += loc_f(1);
        dEdq(3*(n0-1)+2) += loc_f(2);

        // bending force for local node1
        dEdq(3*(n1-1))   += loc_f(3);
        dEdq(3*(n1-1)+1) += loc_f(4);
        dEdq(3*(n1-1)+2) += loc_f(5);

        // bending force for local node2
        dEdq(3*(n2-1))   += loc_f(6);
        dEdq(3*(n2-1)+1) += loc_f(7);
        dEdq(3*(n2-1)+2) += loc_f(8);

        // bending force for local node3
        dEdq(3*(n3-1))   += loc_f(9);
        dEdq(3*(n3-1)+1) += loc_f(10);
        dEdq(3*(n3-1)+2) += loc_f(11);

    }
}

void SolverImpl::calculate_jstretch(Element& el, Eigen::MatrixXd& ddEddq) {

    // loop over 3 edges
    for (int k = 1; k <= 3; k++) {

        double l0 = el.get_len_edge(k);
        int loc1 = 0, loc2 = 0;

        switch (k) {
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

        double ks = m_SimPar->kstretch();

        // local jacobian matrix
        Eigen::MatrixXd loc_j = Eigen::MatrixXd::Zero(6, 6);

        loc_j(0, 0) = ddEs_dx1dx1(loc_xyz1, loc_xyz2, l0, ks);
        loc_j(0, 1) = ddEs_dx1dy1(loc_xyz1, loc_xyz2, l0, ks);
        loc_j(0, 2) = ddEs_dx1dz1(loc_xyz1, loc_xyz2, l0, ks);
        loc_j(0, 3) = ddEs_dx1dx2(loc_xyz1, loc_xyz2, l0, ks);
        loc_j(0, 4) = ddEs_dx1dy2(loc_xyz1, loc_xyz2, l0, ks);
        loc_j(0, 5) = ddEs_dx1dz2(loc_xyz1, loc_xyz2, l0, ks);

        loc_j(1, 1) = ddEs_dy1dy1(loc_xyz1, loc_xyz2, l0, ks);
        loc_j(1, 2) = ddEs_dy1dz1(loc_xyz1, loc_xyz2, l0, ks);
        loc_j(1, 3) = ddEs_dx2dy1(loc_xyz1, loc_xyz2, l0, ks); //
        loc_j(1, 4) = ddEs_dy1dy2(loc_xyz1, loc_xyz2, l0, ks);
        loc_j(1, 5) = ddEs_dy1dz2(loc_xyz1, loc_xyz2, l0, ks);

        loc_j(2, 2) = ddEs_dz1dz1(loc_xyz1, loc_xyz2, l0, ks);
        loc_j(2, 3) = ddEs_dx2dz1(loc_xyz1, loc_xyz2, l0, ks); //
        loc_j(2, 4) = ddEs_dy2dz1(loc_xyz1, loc_xyz2, l0, ks); //
        loc_j(2, 5) = ddEs_dz1dz2(loc_xyz1, loc_xyz2, l0, ks);

        loc_j(3, 3) = ddEs_dx2dx2(loc_xyz1, loc_xyz2, l0, ks);
        loc_j(3, 4) = ddEs_dx2dy2(loc_xyz1, loc_xyz2, l0, ks);
        loc_j(3, 5) = ddEs_dx2dz2(loc_xyz1, loc_xyz2, l0, ks);

        loc_j(4, 4) = ddEs_dy2dy2(loc_xyz1, loc_xyz2, l0, ks);
        loc_j(4, 5) = ddEs_dy2dz2(loc_xyz1, loc_xyz2, l0, ks);
        loc_j(5, 5) = ddEs_dz2dz2(loc_xyz1, loc_xyz2, l0, ks);

        for (int i = 0; i < 2 * m_SimGeo->ndof(); i++) {
            for (int j = 0; j < 2 * m_SimGeo->ndof(); j++) {
                if (i > j)
                    loc_j(i,j) = loc_j(j,i);
            }
        }

        // jacobian needs to be divided by 2, because 1/2 wasn't included during derivation
        loc_j = 0.5 * loc_j;

        // if an edge is shared by two elements, local jacobian will be double counted
        if (el.get_nearby_element(k) != 0)
            loc_j = loc_j * 0.5;

        // local node number corresponds to global node number
        unsigned int n1 = el.get_node_num(loc1);
        unsigned int n2 = el.get_node_num(loc2);

        unsigned int nx1 = 3*(n1-1);
        unsigned int nx2 = 3*(n2-1);

        //std::cout << n1 << '\t' << n2 << '\t' << l0 << std::endl;
        //std::cout << loc_j << std::endl;

        // update jacobian
        for (int i = 0; i < m_SimGeo->ndof(); i++) {
            for (int j = 0; j < m_SimGeo->ndof(); j++) {
                ddEddq(nx1+i, nx1+j) += loc_j(i, j);
                ddEddq(nx1+i, nx2+j) += loc_j(i, j+m_SimGeo->ndof());
                ddEddq(nx2+i, nx1+j) += loc_j(i+m_SimGeo->ndof(), j);
                ddEddq(nx2+i, nx2+j) += loc_j(i+m_SimGeo->ndof(), j+m_SimGeo->ndof());
            }
        }
    }
}

void SolverImpl::calculate_jshear(Element& el, Eigen::MatrixXd& ddEddq) {

    // coordinates of local nodes
    Eigen::Vector3d loc_xyz1 = el.get_node(1)->get_xyz();
    Eigen::Vector3d loc_xyz2 = el.get_node(2)->get_xyz();
    Eigen::Vector3d loc_xyz3 = el.get_node(3)->get_xyz();

    // area
    double a0 = el.get_area();

    double ksh = m_SimPar->kshear();

    // local jacobian matrix
    Eigen::MatrixXd loc_j = Eigen::MatrixXd::Zero(9, 9);

    loc_j(0, 0) = ddEsh_dx1dx1(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(0, 1) = ddEsh_dx1dy1(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(0, 2) = ddEsh_dx1dz1(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(0, 3) = ddEsh_dx1dx2(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(0, 4) = ddEsh_dx1dy2(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(0, 5) = ddEsh_dx1dz2(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(0, 6) = ddEsh_dx1dx3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(0, 7) = ddEsh_dx1dy3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(0, 8) = ddEsh_dx1dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);

    loc_j(1, 1) = ddEsh_dy1dy1(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(1, 2) = ddEsh_dy1dz1(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(1, 3) = ddEsh_dx2dy1(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh); //
    loc_j(1, 4) = ddEsh_dy1dy2(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(1, 5) = ddEsh_dy1dz2(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(1, 6) = ddEsh_dx3dy1(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh); //
    loc_j(1, 7) = ddEsh_dy1dy3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(1, 8) = ddEsh_dy1dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);

    loc_j(2, 2) = ddEsh_dz1dz1(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(2, 3) = ddEsh_dx2dz1(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh); //
    loc_j(2, 4) = ddEsh_dy2dz1(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh); //
    loc_j(2, 5) = ddEsh_dz1dz2(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(2, 6) = ddEsh_dx3dz1(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh); //
    loc_j(2, 7) = ddEsh_dy3dz1(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh); //
    loc_j(2, 8) = ddEsh_dz1dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);

    loc_j(3, 3) = ddEsh_dx2dx2(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(3, 4) = ddEsh_dx2dy2(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(3, 5) = ddEsh_dx2dz2(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(3, 6) = ddEsh_dx2dx3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(3, 7) = ddEsh_dx2dy3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(3, 8) = ddEsh_dx2dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);

    loc_j(4, 4) = ddEsh_dy2dy2(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(4, 5) = ddEsh_dy2dz2(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(4, 6) = ddEsh_dx2dx3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(4, 7) = ddEsh_dx2dy3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(4, 8) = ddEsh_dx2dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);

    loc_j(5, 5) = ddEsh_dz2dz2(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(5, 6) = ddEsh_dx3dz2(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh); //
    loc_j(5, 7) = ddEsh_dy3dz2(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh); //
    loc_j(5, 8) = ddEsh_dz2dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);

    loc_j(6, 6) = ddEsh_dx3dx3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(6, 7) = ddEsh_dx3dy3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(6, 8) = ddEsh_dx3dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);

    loc_j(7, 7) = ddEsh_dy3dy3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(7, 8) = ddEsh_dy3dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);
    loc_j(8, 8) = ddEsh_dz3dz3(loc_xyz1, loc_xyz2, loc_xyz3, a0, ksh);

    for (int i = 0; i < 3 * m_SimGeo->ndof(); i++) {
        for (int j = 0; j < 3 * m_SimGeo->ndof(); j++) {
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
    for (int i = 0; i < m_SimGeo->ndof(); i++) {
        for (int j = 0; j < m_SimGeo->ndof(); j++) {

            ddEddq(nx1+i, nx1+j) += loc_j(i, j);
            ddEddq(nx1+i, nx2+j) += loc_j(i, j+m_SimGeo->ndof());
            ddEddq(nx1+i, nx3+j) += loc_j(i, j+m_SimGeo->ndof()*2);

            ddEddq(nx2+i, nx1+j) += loc_j(i+m_SimGeo->ndof(), j);
            ddEddq(nx2+i, nx2+j) += loc_j(i+m_SimGeo->ndof(), j+m_SimGeo->ndof());
            ddEddq(nx2+i, nx3+j) += loc_j(i+m_SimGeo->ndof(), j+m_SimGeo->ndof()*2);

            ddEddq(nx3+i, nx1+j) += loc_j(i+m_SimGeo->ndof()*2, j);
            ddEddq(nx3+i, nx2+j) += loc_j(i+m_SimGeo->ndof()*2, j+m_SimGeo->ndof());
            ddEddq(nx3+i, nx3+j) += loc_j(i+m_SimGeo->ndof()*2, j+m_SimGeo->ndof()*2);
        }
    }
}

void SolverImpl::calculate_jbend(Element& el, Eigen::MatrixXd& ddEddq) {

    // loop over 3 edges
    for (int i = 1; i <= 3; i++) {

        // if the edge is a hinge
        if ( !el.is_hinge(i) )
            continue;

        // get the pointer to that hinge
        Hinge* temp_hinge = el.get_hinge(i);

        // and has NOT been visited yet
        if (temp_hinge->check_visited())
            continue;

        // mark the edge as visited
        temp_hinge->mark_visited();

        // coordinates of local nodes
        Eigen::Vector3d loc_xyz0 = temp_hinge->get_node(0)->get_xyz();
        Eigen::Vector3d loc_xyz1 = temp_hinge->get_node(1)->get_xyz();
        Eigen::Vector3d loc_xyz2 = temp_hinge->get_node(2)->get_xyz();
        Eigen::Vector3d loc_xyz3 = temp_hinge->get_node(3)->get_xyz();

        // original hinge length
        double l0 = el.get_len_edge(i);

        // original sum of the area of two elements
        double a0 = temp_hinge->get_area_sum();

        // original norm(n1)*norm(n2)
        double n1n20 = temp_hinge->get_n1n20();

        // original kappa0
        Eigen::Vector3d kappa0 = temp_hinge->get_kappa0();

        // local jacobian matrix
        Eigen::MatrixXd loc_j = Eigen::MatrixXd::Zero(12, 12);

        loc_j(0, 0) = ddEb_dx0dx0(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(0, 1) = ddEb_dx0dy0(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(0, 2) = ddEb_dx0dz0(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(0, 3) = ddEb_dx0dx1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(0, 4) = ddEb_dx0dy1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(0, 5) = ddEb_dx0dz1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(0, 6) = ddEb_dx0dx2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(0, 7) = ddEb_dx0dy2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(0, 8) = ddEb_dx0dz2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(0, 9) = ddEb_dx0dx3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(0, 10) = ddEb_dx0dy3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(0, 11) = ddEb_dx0dz3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);

        loc_j(1, 1) = ddEb_dy0dy0(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(1, 2) = ddEb_dy0dz0(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(1, 3) = ddEb_dx1dy0(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20); //
        loc_j(1, 4) = ddEb_dy0dy1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(1, 5) = ddEb_dy0dz1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(1, 6) = ddEb_dx2dy0(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20); //
        loc_j(1, 7) = ddEb_dy0dy2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(1, 8) = ddEb_dy0dz2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(1, 9) = ddEb_dx3dy0(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20); //
        loc_j(1, 10) = ddEb_dy0dy3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(1, 11) = ddEb_dy0dz3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);

        loc_j(2, 2) = ddEb_dz0dz0(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(2, 3) = ddEb_dx1dz0(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20); //
        loc_j(2, 4) = ddEb_dy1dz0(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20); //
        loc_j(2, 5) = ddEb_dz0dz1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(2, 6) = ddEb_dx2dz0(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20); //
        loc_j(2, 7) = ddEb_dy2dz0(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20); //
        loc_j(2, 8) = ddEb_dz0dz2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(2, 9) = ddEb_dx3dz0(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20); //
        loc_j(2, 10) = ddEb_dy3dz0(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);//
        loc_j(2, 11) = ddEb_dz0dz3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);

        loc_j(3, 3) = ddEb_dx1dx1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(3, 4) = ddEb_dx1dy1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(3, 5) = ddEb_dx1dz1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(3, 6) = ddEb_dx1dx2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(3, 7) = ddEb_dx1dy2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(3, 8) = ddEb_dx1dz2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(3, 9) = ddEb_dx1dx3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(3, 10) = ddEb_dx1dy3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(3, 11) = ddEb_dx1dz3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);

        loc_j(4, 4) = ddEb_dy1dy1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(4, 5) = ddEb_dy1dz1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(4, 6) = ddEb_dx2dy1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20); //
        loc_j(4, 7) = ddEb_dy1dy2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(4, 8) = ddEb_dy1dz2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(4, 9) = ddEb_dx3dy1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20); //
        loc_j(4, 10) = ddEb_dy1dy3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(4, 11) = ddEb_dy1dz3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);

        loc_j(5, 5) = ddEb_dz1dz1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(5, 6) = ddEb_dx2dz1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20); //
        loc_j(5, 7) = ddEb_dy2dz1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20); //
        loc_j(5, 8) = ddEb_dz1dz2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(5, 9) = ddEb_dx3dz1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20); //
        loc_j(5, 10) = ddEb_dy3dz1(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);//
        loc_j(5, 11) = ddEb_dz1dz3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);

        loc_j(6, 6) = ddEb_dx2dx2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(6, 7) = ddEb_dx2dy2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(6, 8) = ddEb_dx2dz2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(6, 9) = ddEb_dx2dx3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(6, 10) = ddEb_dx2dy3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(6, 11) = ddEb_dx2dz3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);

        loc_j(7, 7) = ddEb_dy2dy2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(7, 8) = ddEb_dy2dz2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(7, 9) = ddEb_dx3dy2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20); //
        loc_j(7, 10) = ddEb_dy2dy3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(7, 11) = ddEb_dy2dz3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);

        loc_j(8, 8) = ddEb_dz2dz2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(8, 9) = ddEb_dx3dz2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20); //
        loc_j(8, 10) = ddEb_dy3dz2(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);//
        loc_j(8, 11) = ddEb_dz2dz3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);

        loc_j(9, 9) = ddEb_dx3dx3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(9, 10) = ddEb_dx3dy3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(9, 11) = ddEb_dx3dz3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);

        loc_j(10, 10) = ddEb_dy3dy3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);
        loc_j(10, 11) = ddEb_dy3dz3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20); //
        loc_j(11, 11) = ddEb_dz3dz3(loc_xyz0, loc_xyz1, loc_xyz2, loc_xyz3, kappa0, n1n20);

        for (int k = 0; k < 4 * m_SimGeo->ndof(); k++) {
            for (int j = 0; j < 4 * m_SimGeo->ndof(); j++) {
                if (k > j)
                    loc_j(k,j) = loc_j(j,k);
            }
        }

        // force needs to be multiplied by the coefficients
        loc_j = loc_j * 0.5 * m_SimPar->kbend() * 3 * pow(l0, 2) / a0;


        //std::cout << "local bending jacobian" << std::endl;
        //std::cout << loc_j << std::endl;

        // local node number corresponds to global node number
        unsigned int n0 = temp_hinge->get_node_num(0);
        unsigned int n1 = temp_hinge->get_node_num(1);
        unsigned int n2 = temp_hinge->get_node_num(2);
        unsigned int n3 = temp_hinge->get_node_num(3);

        unsigned int nx0 = 3*(n0-1);
        unsigned int nx1 = 3*(n1-1);
        unsigned int nx2 = 3*(n2-1);
        unsigned int nx3 = 3*(n3-1);

        // update jacobian
        for (int k = 0; k < m_SimGeo->ndof(); k++) {
            for (int j = 0; j < m_SimGeo->ndof(); j++) {

                ddEddq(nx0+k, nx0+j) += loc_j(k, j);
                ddEddq(nx0+k, nx1+j) += loc_j(k, j+m_SimGeo->ndof());
                ddEddq(nx0+k, nx2+j) += loc_j(k, j+m_SimGeo->ndof()*2);
                ddEddq(nx0+k, nx3+j) += loc_j(k, j+m_SimGeo->ndof()*3);

                ddEddq(nx1+k, nx0+j) += loc_j(k+m_SimGeo->ndof(), j);
                ddEddq(nx1+k, nx1+j) += loc_j(k+m_SimGeo->ndof(), j+m_SimGeo->ndof());
                ddEddq(nx1+k, nx2+j) += loc_j(k+m_SimGeo->ndof(), j+m_SimGeo->ndof()*2);
                ddEddq(nx1+k, nx3+j) += loc_j(k+m_SimGeo->ndof(), j+m_SimGeo->ndof()*3);

                ddEddq(nx2+k, nx0+j) += loc_j(k+m_SimGeo->ndof()*2, j);
                ddEddq(nx2+k, nx1+j) += loc_j(k+m_SimGeo->ndof()*2, j+m_SimGeo->ndof());
                ddEddq(nx2+k, nx2+j) += loc_j(k+m_SimGeo->ndof()*2, j+m_SimGeo->ndof()*2);
                ddEddq(nx2+k, nx3+j) += loc_j(k+m_SimGeo->ndof()*2, j+m_SimGeo->ndof()*3);

                ddEddq(nx3+k, nx0+j) += loc_j(k+m_SimGeo->ndof()*3, j);
                ddEddq(nx3+k, nx1+j) += loc_j(k+m_SimGeo->ndof()*3, j+m_SimGeo->ndof());
                ddEddq(nx3+k, nx2+j) += loc_j(k+m_SimGeo->ndof()*3, j+m_SimGeo->ndof()*2);
                ddEddq(nx3+k, nx3+j) += loc_j(k+m_SimGeo->ndof()*3, j+m_SimGeo->ndof()*3);
            }
        }
    }
}