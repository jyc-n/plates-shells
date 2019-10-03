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
    //*      hanging 1 corner
    //* ----------------------------
    for (int i = 0; i < 3; i++)
        m_dirichletDofs.push_back(i);

    //* ----------------------------
    //*      hanging 2 corners
    //* ----------------------------
    //    for (int i = 0; i < 3; i++) {
    //        m_dirichletDofs.push_back(i);
    //        m_dirichletDofs.push_back(i + 3*(m_SimGeo->num_nodes_len()-1));
    //    }

    //* ----------------------------
    //*      draping middle
    //* ----------------------------
    //  int mid = m_SimGeo->num_nodes_len() / 2;
    //  int mid_dof = 3 * (m_SimGeo->num_nodes_len() * (mid-1) + mid - 1);
    //  for (int i = 0; i < 18; i++) {
    //      m_dirichletDofs.push_back(i+mid_dof - 6*m_SimGeo->num_nodes_len());
    //      m_dirichletDofs.push_back(i+mid_dof - 3*m_SimGeo->num_nodes_len());
    //      m_dirichletDofs.push_back(i+mid_dof);
    //      m_dirichletDofs.push_back(i+mid_dof + 3*m_SimGeo->num_nodes_len());
    //      m_dirichletDofs.push_back(i+mid_dof + 6*m_SimGeo->num_nodes_len());
    //      m_dirichletDofs.push_back(i+mid_dof + 9*m_SimGeo->num_nodes_len());
    //  }

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