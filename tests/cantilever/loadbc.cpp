#include <iostream>
#include <algorithm>
#include "loadbc.h"
#include "parameters.h"
#include "geometry.h"

// constructor
Boundary::Boundary(Parameters* SimPar, Geometry* SimGeo) {
    m_SimPar = SimPar;
    m_SimGeo = SimGeo;
    ENABLE_GRAVITY = true;
}

// destructor
Boundary::~Boundary() {}

// by default, rectangular plate
//    # of sets for default rectangular plate
//                  4
//           _________________
//        |                   |
//    1   |                   |   3
//        |                   |
//        | _________________ |
//                  2
//
//*     Several test cases has been implementmented here
//*     Uncomment the code to run different test cases
//
void Boundary::initBC() {
    //* ----------------------------
    //*        Force vector
    //* ----------------------------
    m_fext = Eigen::VectorXd::Zero(m_SimGeo->nn() * m_SimGeo->nsd());

    //* ----------------------------
    //*     Gravity option
    //* ----------------------------
    // gravity
    if (ENABLE_GRAVITY) {
        configGravity(Z_DIR, false);
    }

    //* ----------------------------
    //*      cantilever plate
    //* ----------------------------
    //! Note: when using clamped bc, need to modify the number of nodes!
    // clamped edge 1
    Sets edge1(false);
    edge1.setFixed(true, true, true);
    for (int num = 1; num < m_SimGeo->nn(); num += m_SimGeo->num_nodes_len()) {
        edge1.m_nodes.push_back(num);
        edge1.m_nodes.push_back(num+1);
    }
    findDirichletDofs(edge1);
    //* ----------------------------

    //* sort the array for Dirichlet dofs
    std::sort(m_dirichletDofs.begin(), m_dirichletDofs.end());
}

// configure external force
void Boundary::configForce(const Sets& t_set) {
    for (auto inode = t_set.m_nodes.begin(); inode != t_set.m_nodes.end(); inode++) {
        for (int dir = 0; dir <= 2; dir++) {
            m_fext(m_SimGeo->nsd()*(*inode-1) + dir) = t_set.m_force[dir];
        }
    }
}

// enable gravity
void Boundary::configGravity(const int dir, bool sign) {
    double sign_g = (sign) ? 1.0 : -1.0;
    for (int i = dir; i < m_SimGeo->nn() * m_SimGeo->nsd(); i+=3) {
        m_fext(i) = sign_g * m_SimGeo->m_mass(i) * m_SimPar->gconst();
    }
}

void Boundary::findDirichletDofs(const Sets& t_set) {
    for (auto inode = t_set.m_nodes.begin(); inode != t_set.m_nodes.end(); inode++) {
        for (int dir = 0; dir <= 2; dir++) {
            if (t_set.m_fixed[dir]) {
                m_dirichletDofs.push_back(m_SimGeo->nsd()*(*inode-1) + dir);
            }
        }
    }
}

bool Boundary::inDirichletBC(int index) {
    auto iter = find(m_dirichletDofs.begin(), m_dirichletDofs.end(), index);
    return !(iter == m_dirichletDofs.end());
}