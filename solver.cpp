#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

#include "solver.h"
#include "parameters.h"
#include "geometry.h"
#include "loadbc.h"
#include "node.h"
#include "hinge.h"
#include "edge.h"
#include "element.h"
#include "utilities.h"
#include "stretching.h"
#include "shearing.h"
#include "bending.h"


SolverImpl::SolverImpl(Parameters* SimPar, Geometry* SimGeo, Boundary* SimBC) {
    m_SimPar = SimPar;
    m_SimGeo = SimGeo;
    m_SimBC  = SimBC;
}

void SolverImpl::initSolver() {
    m_dEdq = Eigen::VectorXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof());
    m_ddEddq = Eigen::MatrixXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof(), m_SimGeo->nn() * m_SimGeo->ndof());
    m_residual = Eigen::VectorXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof());
    m_jacobian = Eigen::MatrixXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof(), m_SimGeo->nn() * m_SimGeo->ndof());
}

// ========================================= //
//            Main solver function           //
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
    double tol = m_SimGeo->m_mi * 9.81 * 1e-1;

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

            initSolver();

            std::cout << niter+1 << '\t';

            t.start();

            // calculate derivatives of energy functions
            calcDEnergy(m_dEdq, m_ddEddq);

            // calculate residual vector
            updateResidual(dof, dof_new, vel, m_dEdq, m_residual);
            //std::cout << m_residual << std::endl;

            // calculate jacobian matrix
            updateJacobian(m_ddEddq, m_jacobian);
            //std::cout << m_jacobian << std::endl;

            // residual vector and jacobian matrix at dofs that are not specified
            Eigen::VectorXd res_free = unconsVec(m_residual);
            Eigen::MatrixXd jacobian_free = unconsMat(m_jacobian);
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

void SolverImpl::calcDEnergy(Eigen::VectorXd& dEdq, Eigen::MatrixXd& ddEddq) {

    Eigen::VectorXd fstretch = Eigen::VectorXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof());
    Eigen::VectorXd fshear   = Eigen::VectorXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof());
    Eigen::VectorXd fbend    = Eigen::VectorXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof());

    Eigen::MatrixXd jstretch = Eigen::MatrixXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof(), m_SimGeo->nn() * m_SimGeo->ndof());
    Eigen::MatrixXd jshear   = Eigen::MatrixXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof(), m_SimGeo->nn() * m_SimGeo->ndof());
    Eigen::MatrixXd jbend    = Eigen::MatrixXd::Zero(m_SimGeo->nn() * m_SimGeo->ndof(), m_SimGeo->nn() * m_SimGeo->ndof());


    //Timer ts;
    calcStretch(fstretch, jstretch);
    calcShear(fshear, jshear);
    calcBend(fbend, jbend);

    //std::cout << ts.elapsed() << " ms \t";

    dEdq = fstretch + fshear + fbend;
    ddEddq = jstretch + jshear + jbend;

    /*
    std::cout << "stretch" << std::endl;
    std::cout << fstretch << std::endl;
    std::cout << "shear" << std::endl;
    std::cout << fshear << std::endl;
    std::cout << "bend" << std::endl;
    std::cout << fbend << std::endl;
     */

    /*
    std::cout << "stretch" << std::endl;
    std::cout << jstretch << std::endl;
    std::cout << "shear" << std::endl;
    std::cout << jshear << std::endl;
    std::cout << "bend" << std::endl;
    std::cout << jbend << std::endl;
    */
}

/*
 *    f_i = m_i * (q_i(t_n+1) - q_i(t_n)) / dt^2 - m_i * v(t_n) / dt + dE/dq - F_ext
 */
void SolverImpl::updateResidual(Eigen::VectorXd& qn, Eigen::VectorXd& qnew, Eigen::VectorXd& vel, Eigen::VectorXd& dEdq, Eigen::VectorXd& res_f) {

    double dt = m_SimPar->dt();

    for (unsigned int i = 0; i < m_SimGeo->nn() * m_SimGeo->ndof(); i++)
        res_f(i) = m_SimGeo->m_mass(i) * (qnew(i) - qn(i)) / pow(dt,2) - m_SimGeo->m_mass(i) * vel(i)/dt
                 + dEdq(i) - m_SimBC->m_fext(i);
    /*
    std::cout << "dEdq" << std::endl;
    std::cout << dEdq << std::endl;
    std::cout << "residual" << std::endl;
    std::cout << res_f << std::endl;
    */
}

/*
 *    J_ij = m_i / dt^2 * delta_ij + d^2 E / dq_i dq_j
 */
void SolverImpl::updateJacobian(Eigen::MatrixXd& ddEddq, Eigen::MatrixXd& mat_j) {

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

    //denseSolver(temp_j, temp_f, dq_free);
    sparseSolver(temp_j, temp_f, dq_free);
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
//              Solver functions             //
// ========================================= //

void SolverImpl::denseSolver(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x) {
    x = A.llt().solve(b);
}

void SolverImpl::sparseSolver(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x) {
    // definition numerical error
    double epsilon = 1e-16;
    double refVal = 1.0;

    Eigen::SparseMatrix<double> As;
    As = A.sparseView(refVal, epsilon);

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper> CGsolver;
    CGsolver.compute(As);
    if (CGsolver.info() != Eigen::Success)
        throw "decomposition failed";

    x = CGsolver.solve(b);
    if (CGsolver.info() != Eigen::Success)
        throw "solving failed";
}

// ========================================= //
//          Other helper functions           //
// ========================================= //

Eigen::VectorXd SolverImpl::unconsVec(const Eigen::VectorXd& vec) {

    Eigen::VectorXd vec_free = Eigen::VectorXd::Zero(m_SimBC->m_numFree);

    std::vector<int>::const_iterator p = m_SimBC->m_specifiedDof.begin();
    for (int i = 0, j = 0; i < m_SimBC->m_numTotal && p != m_SimBC->m_specifiedDof.end(); i++) {
        if (i == *p) {
            if (p != m_SimBC->m_specifiedDof.end()-1)
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

// stretch energy for each element
void SolverImpl::calcStretch(Eigen::VectorXd& dEdq, Eigen::MatrixXd& ddEddq){

    // loop over the edge list
    for (std::vector<Edge*>::iterator iedge = m_SimGeo->m_edgeList.begin(); iedge != m_SimGeo->m_edgeList.end(); iedge++) {

        // local stretch force: fx1, fy1, fz1, fx2, fy2, fz2
        Eigen::VectorXd loc_f = Eigen::VectorXd::Zero(6);

        // local jacobian matrix
        Eigen::MatrixXd loc_j = Eigen::MatrixXd::Zero(6, 6);

        Stretching EStretch(*iedge);

        EStretch.m_ks = m_SimPar->kstretch();

        EStretch.init();
        EStretch.locStretch(loc_f, loc_j);

        // local node number corresponds to global node number
        unsigned int n1 = (*iedge)->get_node_num(1);
        unsigned int n2 = (*iedge)->get_node_num(2);

        //std::cout << n1 << '\t' << n2 << '\n';
        //std::cout << loc_f << std::endl;
        //std::cout << loc_j << std::endl;

        // stretch force for local node1
        dEdq(3*(n1-1))   += loc_f(0);
        dEdq(3*(n1-1)+1) += loc_f(1);
        dEdq(3*(n1-1)+2) += loc_f(2);

        // stretch force for local node2
        dEdq(3*(n2-1))   += loc_f(3);
        dEdq(3*(n2-1)+1) += loc_f(4);
        dEdq(3*(n2-1)+2) += loc_f(5);

        // update jacobian
        unsigned int nx1 = 3*(n1-1);
        unsigned int nx2 = 3*(n2-1);

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

// shear energy for each element
void SolverImpl::calcShear(Eigen::VectorXd &dEdq, Eigen::MatrixXd &ddEddq){

    for (std::vector<Element>::iterator iel = m_SimGeo->m_elementList.begin(); iel != m_SimGeo->m_elementList.end(); iel++) {

        // local shearing force: fx1, fy1, fz1, fx2, fy2, fz2, fx3, fy3, fz3
        Eigen::VectorXd loc_f = Eigen::VectorXd::Zero(9);

        // local jacobian matrix
        Eigen::MatrixXd loc_j = Eigen::MatrixXd::Zero(9, 9);

        Shearing EShear(&(*iel));

        EShear.m_coeff = 0.5 * m_SimPar->kshear() * (*iel).get_area();

        EShear.init();
        EShear.locShear(loc_f, loc_j);

        // local node number corresponds to global node number
        unsigned int n1 = (*iel).get_node_num(1);
        unsigned int n2 = (*iel).get_node_num(2);
        unsigned int n3 = (*iel).get_node_num(3);

        //std::cout << n1 << '\t' << n2 << '\t' << n3 << '\n';
        //std::cout << loc_f << std::endl;
        //std::cout << loc_j << std::endl;

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

        // update jacobian
        unsigned int nx1 = 3*(n1-1);
        unsigned int nx2 = 3*(n2-1);
        unsigned int nx3 = 3*(n3-1);

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
}

void SolverImpl::calcBend(Eigen::VectorXd &dEdq, Eigen::MatrixXd &ddEddq){

    // loop over the hinge list
    for (std::vector<Hinge*>::iterator ihinge = m_SimGeo->m_hingeList.begin(); ihinge != m_SimGeo->m_hingeList.end(); ihinge++) {

        // local bending force: fx0, fy0, fz0, fx1, fy1, fz1, fx2, fy2, fz2, fx3, fy3, fz3
        Eigen::VectorXd loc_f = Eigen::VectorXd::Zero(12);

        // local jacobian matrix
        Eigen::MatrixXd loc_j = Eigen::MatrixXd::Zero(12, 12);

        // bending energy calculation
        Bending Ebend(*ihinge);

        Ebend.m_coeff = (*ihinge)->m_const * m_SimPar->kbend();

        Ebend.init();
        Ebend.locBend(loc_f, loc_j);

        // local node number corresponds to global node number
        unsigned int n0 = (*ihinge)->get_node_num(0);
        unsigned int n1 = (*ihinge)->get_node_num(1);
        unsigned int n2 = (*ihinge)->get_node_num(2);
        unsigned int n3 = (*ihinge)->get_node_num(3);

/*
        std::cout << n0 << "\t" << n1 << "\t" << n2 << "\t" << n3 << std::endl;

        std::cout << "m_k " << (*ihinge)->m_k << std::endl;
        std::cout << loc_j << std::endl;*/


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

        // update jacobian
        unsigned int nx0 = 3*(n0-1);
        unsigned int nx1 = 3*(n1-1);
        unsigned int nx2 = 3*(n2-1);
        unsigned int nx3 = 3*(n3-1);

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