#include <iostream>
#include <iomanip>
#include <fstream>

#include "geometry.h"
#include "parameters.h"
#include "node.h"
#include "element.h"
#include "edge.h"
#include "hinge.h"

#include <Eigen/Dense>

Geometry::Geometry(Parameters* SimPar) {
    m_SimPar  = SimPar;
}

Geometry::~Geometry() {
    for (std::vector<Edge*>::iterator iedge = m_edgeList.begin(); iedge != m_edgeList.end(); iedge++)
        delete *iedge;
    for (std::vector<Hinge*>::iterator ihinge = m_hingeList.begin(); ihinge != m_hingeList.end(); ihinge++)
        delete *ihinge;
}

void Geometry::buildGeo() {
    m_coord       = Eigen::MatrixXd::Zero(m_nn, m_ndof);
    m_conn        = Eigen::MatrixXi::Zero(m_nel, m_nen);
    m_dof         = Eigen::VectorXd::Zero(m_nn * m_ndof);
    map_nodes     = Eigen::MatrixXi::Zero(m_num_nodes_wid, m_num_nodes_len);
}

void Geometry::calcMass() {
    m_mass = Eigen::VectorXd::Zero(m_nn * m_ndof);

    // Total mass per small rectangle (2 triangular elements)
    double mi = (m_rec_len * m_rec_wid * m_SimPar->thk() * m_SimPar->rho())
                / ((m_num_nodes_len - 1) * (m_num_nodes_wid - 1));

    for (int i = 0; i < m_num_nodes_len; i++) {
        for (int j = 0; j < m_num_nodes_wid; ++j) {

            // position in the mass vector
            int k = i + j * m_num_nodes_len;

            // corner nodes
            if ( (i == 0 || i == m_num_nodes_len - 1) && (j == 0 || j == m_num_nodes_wid - 1) ) {
                m_mass(k * m_ndof) = mi / 4.0;
                m_mass(k * m_ndof + 1) = mi / 4.0;
                m_mass(k * m_ndof + 2) = mi / 4.0;
            }
            // edge nodes
            else if ( (i == 0 || i == m_num_nodes_len - 1) || (j == 0 || j == m_num_nodes_wid - 1) ) {
                m_mass(k * m_ndof) = mi / 2.0;
                m_mass(k * m_ndof + 1) = mi / 2.0;
                m_mass(k * m_ndof + 2) = mi / 2.0;
            }
            // middle nodes
            else {
                m_mass(k * m_ndof) = mi;
                m_mass(k * m_ndof + 1) = mi;
                m_mass(k * m_ndof + 2) = mi;
            }
        }
    }
}

// print geometric information and generate connectivity file
void Geometry::printGeo() {
    std::cout << "----List of Nodes----" << '\n';
    for (int i = 0; i < m_nn; i++) {
        std::cout << i+1 << '\t';
        for (int j = 0; j < m_ndof; j++) {
            std::cout << std::setprecision(8) << std::fixed << m_coord(i,j) << '\t';
        }
        std::cout << std::endl;
    }
    //std::cout << "----DOF Vector----" << '\n';
    //std::cout << m_dof << std::endl;
    //std::cout << "----Map of Nodes----" << '\n';
    //std::cout << map_nodes << std::endl;
    std::cout << "----List of Elements----" << '\n';
    std::string filepath = m_SimPar->outputPath();
    std::string filename = "connectivity.txt";
    std::ofstream myfile((filepath+filename).c_str());

    for (int i = 0; i < m_nel; i++) {
        std::cout << i+1 << '\t';
        myfile << i+1 << '\t';
        for (int j = 0; j < m_nen; j++) {
            std::cout << m_conn(i,j) << '\t';
            myfile << m_conn(i,j) << '\t';
        }
        std::cout << std::endl;
        myfile << std::endl;
    }
}