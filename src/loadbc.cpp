#include <iostream>
#include "loadbc.h"
#include "parameters.h"
#include "geometry.h"

// constructor
Boundary::Boundary(Parameters* SimPar, Geometry* SimGeo) {
    m_SimPar = SimPar;
    m_SimGeo = SimGeo;
    m_2D = false;
}

// destructor
Boundary::~Boundary() {
    for (int i = 0; i < m_sets.size(); i++)
        delete m_sets[i];
}


void Boundary::initBC() {
    getSets();
    configBC();
    //buildBCinfo();
    findDirichletDofs();
}

// by default, rectangular plate
/*    # of sets for default rectangular plate
 *                    4
 *             _________________
 *          |                   |
 *      1   |                   |   3
 *          |                   |
 *          | _________________ |
 *                    2
 */
void Boundary::getSets() {
    m_bounds = 4;

    for (int i = 0; i < m_bounds; i++)
        m_sets.push_back(new Sets);

    for (int num = 1; num < m_SimGeo->nn(); num += m_SimGeo->num_nodes_len())
        m_sets[0]->m_nodes.push_back(num);

    for (int num = 1; num <= m_SimGeo->num_nodes_len(); num++ )
        m_sets[1]->m_nodes.push_back(num);

    for (int num = m_SimGeo->num_nodes_len(); num <= m_SimGeo->nn(); num += m_SimGeo->num_nodes_len())
        m_sets[2]->m_nodes.push_back(num);

    for (int num = m_sets[0]->m_nodes[m_SimGeo->num_nodes_wid()-1]; num <= m_SimGeo->nn(); num++)
        m_sets[3]->m_nodes.push_back(num);
}

// configure the boundary conditions
void Boundary::configBC() {

    m_disp.resize(m_SimGeo->nn() * m_SimGeo->nsd(), true);
    m_fext = Eigen::VectorXd::Zero(m_SimGeo->nn() * m_SimGeo->nsd());

    /*
     *  Uniaxial loading
     *
    m_2D = true;
    m_sets[0]->m_free = false;          // pin left edge
    m_sets[0]->m_typeBC = 1;
    m_sets[2]->m_free = false;          // apply force on the right edge
    m_sets[2]->m_force[0] = 100;
    m_disp[1] = false;                  // fix y dof of node 1
     */

    /*
     *   Simple bending test
     *
    m_sets[0]->m_free = false;          // clamp left edge
    m_sets[0]->m_typeBC = 3;
    fixDir(2);                          // fix y dof of all nodes
    m_sets[2]->m_free = false;          // apply force on the right edge
    m_sets[2]->m_force[2] = -10;
    */

    /*
     *   Bending under gravity
     *
     */
    m_sets[0]->m_free = false;          // clamp left edge
    m_sets[0]->m_typeBC = 3;

    for (int i = 0; i < m_SimGeo->nn() * m_SimGeo->nsd(); i++) {
        if (i % 3 == 2)
            m_fext(i) = - m_SimGeo->m_mass(i) * m_SimPar->gconst();
    }

}

void Boundary::fixDir(int num) {
    switch (num) {
        case 1:
            for (int i = 1; i <= m_SimGeo->nn(); i++)
                m_disp[3*(i-1)] = false;
            break;
        case 2:
            for (int i = 1; i <= m_SimGeo->nn(); i++)
                m_disp[3*(i-1)+1] = false;
            break;
        case 3:
            for (int i = 1; i <= m_SimGeo->nn(); i++)
                m_disp[3*(i-1)+2] = false;
            break;
        default:
            throw "wrong direction!";
    }
}

// build force vector and find specified dofs for the solver
void Boundary::buildBCinfo() {

    if (m_2D)                          // 2D case
        fixDir(3);


    for (const auto& i : m_sets) {

        if ( i->m_free )               // if NO BC applied to this set
            continue;

        if ( i->m_typeBC != 0 ) {       // displacement BC
            switch (i->m_typeBC) {
                case 1:                 // pinned in x
                    for (const auto& nn : i->m_nodes)
                        m_disp[3*(nn-1)] = false;
                    break;
                case 2:                 // pinned in y
                    for (const auto& nn : i->m_nodes)
                        m_disp[3*(nn-1)+1] = false;
                    break;
                case 3:                 // clamped
                    for (const auto& nn : i->m_nodes) {
                        m_disp[3*(nn-1)] = false;
                        m_disp[3*(nn-1)+1] = false;
                        m_disp[3*(nn-1)+2] = false;
                    }

                    // TODO: implement a better algorithm to configure clamped BC
                    for (int num = 1+m_SimGeo->num_nodes_len(); num <= 2*m_SimGeo->num_nodes_len(); num++) {
                        m_disp[3*(num-1)] = false;
                        m_disp[3*(num-1)+1] = false;
                        m_disp[3*(num-1)+2] = false;
                    }
                    break;
                case 4:
                    throw "Not implemented yet";
            }
        }
        else {                          // force BC
            for (const auto& nn : i->m_nodes)
                for (int k = 0; k < 3; k++)
                    m_fext(3*(nn-1)+k) = i->m_force[k];
        }
    }
}

void Boundary::findDirichletDofs() {


    // Hanging cloth corner
    for (int i = 0; i < 6; i++)
        m_dirichletDofs.push_back(i);
    for (int i = m_SimGeo->num_nodes_len()*3; i < m_SimGeo->num_nodes_len()*3+6; i++)
        m_dirichletDofs.push_back(i);
    /*
    // clamped bottom edge
    for (int i = 1; i <= m_SimGeo->num_nodes_len() * 2; i++) {
        for (int j = 0; j < 3; j++)
            m_dirichletDofs.push_back(3*(i-1)+j);
    }

    // Hanging cloth middle
        int mid = m_SimGeo->num_nodes_len() / 2;
        for (int i = mid; i < mid+6; i++)
            m_dirichletDofs.push_back(i);
    */
    //for (auto p = m_dirichletDofs.begin(); p != m_dirichletDofs.end(); p++)
    //    std::cout << *p << '\n';
}

bool Boundary::inDirichletBC(int index) {
    auto iter = find(m_dirichletDofs.begin(), m_dirichletDofs.end(), index);
    return !(iter == m_dirichletDofs.end());
}