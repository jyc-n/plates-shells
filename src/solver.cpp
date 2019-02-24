#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cmath>

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
    DYNAMIC_SOLVER = m_SimPar->solver_op();
    WRITE_OUTPUT = m_SimPar->outop();

    m_numTotal = m_SimGeo->nn() * m_SimGeo->nsd();
    m_numDirichlet = (int) m_SimBC->m_dirichletDofs.size();
    m_numNeumann = m_numTotal - m_numDirichlet;
    
    m_tol = m_SimGeo->m_mi * m_SimPar->gconst() * m_SimPar->ctol();
    m_incRatio = 1.0 / ((double) m_SimPar->nst());

    findMappingVectors();
}

// ========================================= //
//            Main solver function           //
// ========================================= //

// TODO: implement Newmark-beta method
// Time stepping using backward Euler
// FIXME: should consider the acceleration!
bool SolverImpl::step(const int ist, VectorNodes& x, VectorNodes& x_new, VectorNodes& vel) {
    std::cout << "--------Step " << ist << "--------" << std::endl;

    // apply Newton-Raphson Method
    for (int niter = 0; niter < m_SimPar->iter_lim(); niter++) {
        Timer t;

        VectorN rhs(m_numNeumann); rhs.fill(0.0);
        SparseEntries entries_full;

        //Timer t1;
        // calculate derivatives of energy functions
        VectorN dEdq(m_numTotal); dEdq.fill(0.0);
        findDEnergy(dEdq, entries_full);
        //std::cout << t1.elapsed() << '\t';

        // calculate residual vector
        findResidual(vel, x, x_new, dEdq, rhs);

        // display residual
        double error = rhs.norm();
        std::cout << "iter" << niter+1 << '\t' << "error = " << error << '\t';

        // check convergence
        if (error < m_tol) {
            // calculate new velocity vector
            auto findVel = [&x, &x_new] (double dt, VectorNodes& vel) {
                for (int i = 0; i < x.size(); i++)
                    vel[i] = (x_new[i] - x[i]) / dt;
            };
            findVel(m_SimPar->dt(), vel);
            // save current nodal position for next step
            for (int i = 0; i < m_SimGeo->nn(); i++)
                x[i] = x_new[i];
            // display iteration time
            std::cout << "t_iter = " << t.elapsed() << " ms" << std::endl;
            return true;
        }

        // calculate jacobian matrix
        SpMatrix jacobian(m_numNeumann, m_numNeumann);
        findJacobian(entries_full, jacobian);

        // solve for new dof vector
        findDofnew(rhs, jacobian, x_new);

        // display iteration time
        std::cout << "t_iter = " << t.elapsed() << " ms" << std::endl;
    }
    return false;
}

void SolverImpl::dynamic() {
    if (SOLVER_TYPE == 0)
        std::cout << "Eigen CG solver will be used" << std::endl;
    else if (SOLVER_TYPE == 1)
        std::cout << "Pardiso solver will be used" << std::endl;
    
    Timer t_all(true);

    // vector of nodal position t_{n+1}
    VectorNodes nodes_curr(m_SimGeo->nn());
    assert(nodes_curr.size() == m_SimGeo->m_nodes.size());

    // initial guess
    for (int i = 0; i < m_SimGeo->nn(); i++)
        nodes_curr[i] = m_SimGeo->m_nodes[i];

    //vector of nodal velocity
    VectorNodes vel(m_SimGeo->nn());
    for (int i = 0; i < m_SimGeo->nn(); i++)
        vel[i].fill(0.0);
    
    for (int ist = 1; ist <= m_SimPar->nst(); ist++) {
        // stepping
        if (!step(ist, nodes_curr, m_SimGeo->m_nodes, vel)) {
            std::cerr << "Solver did not converge in " << m_SimPar->iter_lim()
                        << " iterations at step " << ist << std::endl;
            throw "Cannot converge! Program terminated";
        }
        if (WRITE_OUTPUT)
            writeToFiles(ist);
    }
    std::cout << "---------------------------" << std::endl;
    std::cout << "Simulation completed" << std::endl;
    std::cout << "Total time used " << t_all.elapsed(true) << " seconds" <<  std::endl;
}

void SolverImpl::quasistatic() {
    // TODO: reimplement quasistatic
}

// write data to files
void SolverImpl::writeToFiles(const int ist) {
    if (ist % m_SimPar->out_freq() != 0)
        return;
    if (m_SimPar->out_freq() == -1)
        if (ist != m_SimPar->nst())
            return;

    // output files
    std::string filepath = m_SimPar->outputPath();
    std::string filename;
    char buffer[20] = {0};
    sprintf(buffer, "result%05d.txt", ist);
    filename.assign(buffer);
    std::ofstream myfile((filepath+filename).c_str());
    for (int k = 0; k < m_SimGeo->nn(); k++) {
        myfile << std::setprecision(8) << std::fixed
                << m_SimGeo->m_nodes[k][0] << '\t'
                << m_SimGeo->m_nodes[k][1] << '\t'
                << m_SimGeo->m_nodes[k][2] << std::endl;
    }
}

// ========================================= //
//       Implementation of subroutines       //
// ========================================= //

void SolverImpl::findDEnergy(VectorN& dEdq, SparseEntries& entries_full) {
    DEStretch(dEdq, entries_full);
    DEShear(dEdq, entries_full);
    DEBend(dEdq, entries_full);
}

/*   TODO: static version will be replaced
 *    f_i = dE/dq - F_ext

void SolverImpl::updateResidual(Eigen::VectorXd& qn, Eigen::VectorXd& qnew, double ratio, Eigen::VectorXd& dEdq, Eigen::VectorXd& res_f) {

    double dt = m_SimPar->dt();

    for (unsigned int i = 0; i < m_SimGeo->nn() * m_SimGeo->nsd(); i++)
        res_f(i) = dEdq(i) - m_SimBC->m_fext(i) * ratio;
}
 */

//
//  f_i = m_i * (q_i(t_n+1) - q_i(t_n)) / dt^2 - m_i * v(t_n) / dt + dE/dq - F_ext
//
void SolverImpl::findResidual(const VectorNodes& vel, const VectorNodes& x, const VectorNodes& x_new, const VectorN& dEdq, VectorN& rhs) {
    double dt = m_SimPar->dt();

    // TODO: better model of viscous damping
    double nu = m_SimPar->vis();
    double area = 0.5 * m_SimGeo->rec_len() * m_SimGeo->rec_wid()
                  / (m_SimGeo->num_nodes_len() * m_SimGeo->num_nodes_wid());

    auto vis = [&dt, &nu, &area] (double qnew, double qn) {
        return nu * area* (qnew - qn) / dt;
    };

    // only take the entries that are NOT in Dirichlet BC
    // CAUTION: if Dirichlet BC is nonzero, need to consider the influence of the Dirichlet BC
    // FIXME: adding Dirichlet influence
    for (int i = 0; i < m_SimGeo->nn(); i++) {
        for (int j = 0; j < m_SimGeo->nsd(); j++) {
            // corresponding position in the full force vector
            int pos = i * m_SimGeo->nsd() + j;
            // check if this is in Neumann BC
            int pos_dof = m_fullToDofs[pos];
            if (pos_dof != -1) {
                rhs(pos_dof) = m_SimGeo->m_mass(pos) * (x_new[i][j] - x[i][j]) / (dt*dt) 
                               - m_SimGeo->m_mass(pos) * vel[i][j]/dt + vis(x_new[i][j], x[i][j])
                               + dEdq(pos) - m_SimBC->m_fext(pos);
            }
        }
    }
}

//
//  J_ij = m_i / dt^2 * delta_ij + d^2 E / dq_i dq_j
//
void SolverImpl::findJacobian(SparseEntries& entries_full, SpMatrix& jacobian) {
    if (DYNAMIC_SOLVER) {
        double dt = m_SimPar->dt();
        // TODO: better model of viscous damping
        double nu = m_SimPar->vis();
        double area = 0.5 * m_SimGeo->rec_len() * m_SimGeo->rec_wid()
                      / (m_SimGeo->num_nodes_len() * m_SimGeo->num_nodes_wid());
        for (int i = 0; i < m_numTotal; i++) {
            // inertia term
            double inertia = m_SimGeo->m_mass(i) / (dt*dt);
            // viscous term
            double viscous = nu * area / dt;
            entries_full.emplace_back(Eigen::Triplet<double>(i, i, inertia + viscous));
        }
    }
    SparseEntries entries_dof;
    for (int p = 0; p < entries_full.size(); p++) {
        int ifull = entries_full[p].row(), jfull = entries_full[p].col();
        int idof = m_fullToDofs[ifull], jdof = m_fullToDofs[jfull];
        if (idof != -1 && jdof != -1) {
            // NOTE: 
            // for Pardiso solver, only store the upper triangular part of jacobian!!
            if (SOLVER_TYPE == 1) {
                // skip the lower triangular part
                if (idof > jdof)
                    continue;
            }
            entries_dof.emplace_back(Eigen::Triplet<double>(idof, jdof, entries_full[p].value()));
        }
    }
    jacobian.setFromTriplets(entries_dof.begin(), entries_dof.end());
}

//
//  q_{n+1} = q_n - J \ f
//
void SolverImpl::findDofnew(const VectorN& rhs, SpMatrix& jacobian, VectorNodes& x_new) {
    VectorN dq(m_numNeumann); dq.fill(0.0);

    //Timer t1;
    if (SOLVER_TYPE == 0)
        sparseSolver(jacobian, rhs, dq);
    else if (SOLVER_TYPE == 1)
        pardisoInterface(jacobian, rhs, dq);
    //std::cout << t1.elapsed() << '\t';

    // map the free part dof vector back to the full dof vector
    // CAUTION: if Dirichlet BC is nonzero, need to consider the motion of the Dirichlet BC
    // FIXME: adding Dirichlet part
    for (int i_dq = 0; i_dq < m_numNeumann; i_dq++) {
        int iN = m_dofsToFull[i_dq] / m_SimGeo->nsd();
        int jN = m_dofsToFull[i_dq] - iN * m_SimGeo->nsd();
        x_new[iN][jN] -= dq(i_dq);
    }
}

// ========================================= //
//              Solver functions             //
// ========================================= //

// conjugate gradient solver from Eigen
void SolverImpl::sparseSolver(const SpMatrix& A, const VectorN& rhs, VectorN& u) {
    Eigen::ConjugateGradient<SpMatrix, Eigen::Upper> CGsolver;
    CGsolver.compute(A);
    if (CGsolver.info() != Eigen::Success)
        throw "decomposition failed";

    u = CGsolver.solve(rhs);
    if (CGsolver.info() != Eigen::Success)
        throw "solving failed";
}

// Pardiso solver interface
void SolverImpl::pardisoInterface(SpMatrix& A, const VectorN& rhs, VectorN& u) {
    // transform the sparse matrix to compressed view
    A.makeCompressed();

    int n = m_numNeumann;
    int nonzeros = A.nonZeros();
    // int innerSize = A.innerSize();
    // int outerSize = A.outerSize();

    /* Matrix data. */
    int* ia = new int[n+1];
    int* ja = new int[nonzeros];
    double* a = new double[nonzeros];

    /* RHS and solution vectors. */
    double* b = new double[n];
    double* x = new double[n];

    // Construct CRS info from Eigen
    // Convert matrix from 0-based C-notation to Fortran 1-based notation
    for (int i = 0; i < n+1; i++)
        ia[i] = A.outerIndexPtr()[i] + 1;
    
    for (int i = 0; i < nonzeros; i++)
        ja[i] = A.innerIndexPtr()[i] + 1;

    for (int i = 0; i < nonzeros; i++)
        a[i] = A.valuePtr()[i];

    for (int i = 0; i < n; i++)
        b[i] = rhs(i);


    // int      nnz = ia[n];
    int      mtype = -2;        /* Real symmetric matrix */
    int      nrhs = 1;          /* Number of right hand sides. */

    /* Internal solver memory pointer pt,                  */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
    /* or void *pt[64] should be OK on both architectures  */ 
    void    *pt[64]; 

    /* Pardiso control parameters. */
    int      iparm[64];
    double   dparm[64];
    int      maxfct, mnum, phase, error, msglvl, solver;

    /* Number of processors. */
    int      num_procs;

    /* Auxiliary variables. */
    char    *var;
    int      i;

    double   ddum;              /* Double dummy */
    int      idum;              /* Integer dummy. */

   
    /* -------------------------------------------------------------------- */
    /* ..  Setup Pardiso control parameters.                                */
    /* -------------------------------------------------------------------- */
    error = 0;
    solver = 0; /* use sparse direct solver */
    pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error); 

    if (error != 0) 
    {
        if (error == -10 )
           printf("No license file found \n");
        if (error == -11 )
           printf("License is expired \n");
        if (error == -12 )
           printf("Wrong username or hostname \n");
        exit(1);
    }
    
    /* Numbers of processors, value of OMP_NUM_THREADS */
    var = getenv("OMP_NUM_THREADS");
    if(var != NULL)
        sscanf( var, "%d", &num_procs );
    else {
        printf("Set environment OMP_NUM_THREADS to 1");
        exit(1);
    }
    iparm[2]  = num_procs;

    maxfct = 1;		/* Maximum number of numerical factorizations.  */
    mnum   = 1;         /* Which factorization to use. */
    
    msglvl = 0;         /* Print statistical information  */
    error  = 0;         /* Initialize error flag */

    //* -------------------------------------------------------------------- */
    /* ..  Reordering and Symbolic Factorization.  This step also allocates */
    /*     all memory that is necessary for the factorization.              */
    /* -------------------------------------------------------------------- */
    phase = 11; 

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
	     &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error, dparm);
    
    if (error != 0) {
        printf("\nERROR during symbolic factorization: %d", error);
        exit(1);
    }
   
    /* -------------------------------------------------------------------- */
    /* ..  Numerical factorization.                                         */
    /* -------------------------------------------------------------------- */    
    phase = 22;
    iparm[32] = 1; /* compute determinant */

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);
   
    if (error != 0) {
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }

    /* -------------------------------------------------------------------- */    
    /* ..  Back substitution and iterative refinement.                      */
    /* -------------------------------------------------------------------- */    
    phase = 33;

    iparm[7] = 1;       /* Max numbers of iterative refinement steps. */
   
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs,
             iparm, &msglvl, b, x, &error,  dparm);
   
    if (error != 0) {
        printf("\nERROR during solution: %d", error);
        exit(3);
    }
    
    for (i = 0; i < n; i++) {
        u(i) = x[i];
    }

    /* -------------------------------------------------------------------- */    
    /* ..  Termination and release of memory.                               */
    /* -------------------------------------------------------------------- */    
    phase = -1;                 /* Release internal memory. */
    
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);

    delete [] ia;
    delete [] ja;
    delete [] a;
    delete [] x;
    delete [] b;
}

// ========================================= //
//          Other helper functions           //
// ========================================= //

// find mapping vectors before the time loop starts
void SolverImpl::findMappingVectors() {
    m_fullToDofs.resize(m_numTotal, -1);
    for (int index = 0; index < m_numTotal; index++) {
        if (!m_SimBC->inDirichletBC(index)) {
            m_dofsToFull.push_back(index);
            m_fullToDofs[index] = (int) (m_dofsToFull.size()-1);
        }
    }
}

// stretch energy for each element
void SolverImpl::DEStretch(VectorN& dEdq, SparseEntries& entries_full) {

    // loop over the edge list
    for (std::vector<Edge*>::iterator iedge = m_SimGeo->m_edgeList.begin(); iedge != m_SimGeo->m_edgeList.end(); iedge++) {

        // local stretch force: fx1, fy1, fz1, fx2, fy2, fz2
        Eigen::VectorXd loc_f = Eigen::VectorXd::Zero(6);
        // local jacobian matrix
        Eigen::MatrixXd loc_j = Eigen::MatrixXd::Zero(6, 6);

        Stretching EStretch(*iedge, m_SimPar->E_modulus(), m_SimPar->thk());
        EStretch.locStretch(loc_f, loc_j);

        // TODO: combine these
        // local node number corresponds to global node number
        unsigned int n1 = (*iedge)->get_node_num(1);
        unsigned int n2 = (*iedge)->get_node_num(2);
        unsigned int nx1 = 3*(n1-1);
        unsigned int nx2 = 3*(n2-1);

        // stretching force
        for (int p = 0; p < 3; p++) {
            dEdq(nx1+p) += loc_f(p);
            dEdq(nx2+p) += loc_f(p+3);
        }

        // stretching jacobian
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                entries_full.emplace_back(Eigen::Triplet<double>(nx1+i, nx1+j, loc_j(i,j)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx1+i, nx2+j, loc_j(i,j+3)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx2+i, nx1+j, loc_j(i+3,j)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx2+i, nx2+j, loc_j(i+3,j+3)));
            }
        }
    }
}

// shear energy for each element
void SolverImpl::DEShear(VectorN& dEdq, SparseEntries& entries_full) {

    // loop over element list
    for (std::vector<Element>::iterator iel = m_SimGeo->m_elementList.begin(); iel != m_SimGeo->m_elementList.end(); iel++) {

        // local shearing force: fx1, fy1, fz1, fx2, fy2, fz2, fx3, fy3, fz3
        Eigen::VectorXd loc_f = Eigen::VectorXd::Zero(9);
        // local jacobian matrix
        Eigen::MatrixXd loc_j = Eigen::MatrixXd::Zero(9, 9);

        Shearing EShear(&(*iel), m_SimPar->E_modulus(), m_SimPar->nu(), (*iel).get_area(), m_SimPar->thk());
        EShear.locShear(loc_f, loc_j);

        // TODO: combine these
        // local node number corresponds to global node number
        unsigned int n1 = (*iel).get_node_num(1);
        unsigned int n2 = (*iel).get_node_num(2);
        unsigned int n3 = (*iel).get_node_num(3);
        unsigned int nx1 = 3*(n1-1);
        unsigned int nx2 = 3*(n2-1);
        unsigned int nx3 = 3*(n3-1);

        // shearing force
        for (int p = 0; p < 3; p++) {
            dEdq(nx1+p) += loc_f(p);
            dEdq(nx2+p) += loc_f(p+3);
            dEdq(nx3+p) += loc_f(p+6);
        }

        // shearing jacobian
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                entries_full.emplace_back(Eigen::Triplet<double>(nx1+i, nx1+j, loc_j(i,j)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx1+i, nx2+j, loc_j(i,j+3)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx1+i, nx3+j, loc_j(i,j+6)));

                entries_full.emplace_back(Eigen::Triplet<double>(nx2+i, nx1+j, loc_j(i+3,j)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx2+i, nx2+j, loc_j(i+3,j+3)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx2+i, nx3+j, loc_j(i+3,j+6)));

                entries_full.emplace_back(Eigen::Triplet<double>(nx3+i, nx1+j, loc_j(i+6,j)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx3+i, nx2+j, loc_j(i+6,j+3)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx3+i, nx3+j, loc_j(i+6,j+6)));
            }
        }
    }
}

void SolverImpl::DEBend(VectorN& dEdq, SparseEntries& entries_full) {

    // loop over the hinge list
    for (std::vector<Hinge*>::iterator ihinge = m_SimGeo->m_hingeList.begin(); ihinge != m_SimGeo->m_hingeList.end(); ihinge++) {

        // local bending force: fx0, fy0, fz0, fx1, fy1, fz1, fx2, fy2, fz2, fx3, fy3, fz3
        Eigen::VectorXd loc_f = Eigen::VectorXd::Zero(12);
        // local jacobian matrix
        Eigen::MatrixXd loc_j = Eigen::MatrixXd::Zero(12, 12);

        // coefficients k for gradient and hessian
        (*ihinge)->m_k = (*ihinge)->m_const * m_SimPar->kbend();

        // bending energy calculation
        Bending Ebend(*ihinge);
        Ebend.locBend(loc_f, loc_j);

        // TODO: combine these
        // local node number corresponds to global node number
        unsigned int n0 = (*ihinge)->get_node_num(0);
        unsigned int n1 = (*ihinge)->get_node_num(1);
        unsigned int n2 = (*ihinge)->get_node_num(2);
        unsigned int n3 = (*ihinge)->get_node_num(3);
        unsigned int nx0 = 3*(n0-1);
        unsigned int nx1 = 3*(n1-1);
        unsigned int nx2 = 3*(n2-1);
        unsigned int nx3 = 3*(n3-1);

        // bending force
        for (int p = 0; p < 3; p++) {
            dEdq(nx0+p) += loc_f(p);
            dEdq(nx1+p) += loc_f(p+3);
            dEdq(nx2+p) += loc_f(p+6);
            dEdq(nx3+p) += loc_f(p+9);
        }

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                entries_full.emplace_back(Eigen::Triplet<double>(nx0+i, nx0+j, loc_j(i,j)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx0+i, nx1+j, loc_j(i,j+3)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx0+i, nx2+j, loc_j(i,j+6)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx0+i, nx3+j, loc_j(i,j+9)));

                entries_full.emplace_back(Eigen::Triplet<double>(nx1+i, nx0+j, loc_j(i+3,j)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx1+i, nx1+j, loc_j(i+3,j+3)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx1+i, nx2+j, loc_j(i+3,j+6)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx1+i, nx3+j, loc_j(i+3,j+9)));

                entries_full.emplace_back(Eigen::Triplet<double>(nx2+i, nx0+j, loc_j(i+6,j)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx2+i, nx1+j, loc_j(i+6,j+3)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx2+i, nx2+j, loc_j(i+6,j+6)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx2+i, nx3+j, loc_j(i+6,j+9)));

                entries_full.emplace_back(Eigen::Triplet<double>(nx3+i, nx0+j, loc_j(i+9,j)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx3+i, nx1+j, loc_j(i+9,j+3)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx3+i, nx2+j, loc_j(i+9,j+6)));
                entries_full.emplace_back(Eigen::Triplet<double>(nx3+i, nx3+j, loc_j(i+9,j+9)));
            }
        }
    }
}

void SolverImpl::analyticalStatic() {
    int nl = m_SimGeo->num_nodes_len();
    double b = m_SimGeo->rec_wid();
    double l = m_SimGeo->rec_len() * double(nl-1)/double(nl);
    double h = m_SimPar->thk();
    double g = 9.81;
    double E = m_SimPar->E_modulus();
    double rho = m_SimPar->rho();

    double mtotal = b*l*h*rho;
    double w = mtotal*g/b;
    double I = 1.0/12.0 * b * pow(h,3);
    double ymax = - w*pow(l,4)/(8.0*E*I);
    std::cout << ymax << std::endl;
}