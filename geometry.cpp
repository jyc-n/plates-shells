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

Geometry::Geometry(Parameters* SimPar)
    : m_datum(0), m_nsd(0), m_nen(0), m_rec_len(0), m_rec_wid(0), m_num_nodes_len(0), m_num_nodes_wid(0),
      m_nn(0), m_nel(0), m_nedge(0), m_nhinge(0)
{
    m_SimPar  = SimPar;
}

Geometry::~Geometry() {
    for (std::vector<Edge*>::iterator iedge = m_edgeList.begin(); iedge != m_edgeList.end(); iedge++)
        delete *iedge;
    for (std::vector<Hinge*>::iterator ihinge = m_hingeList.begin(); ihinge != m_hingeList.end(); ihinge++)
        delete *ihinge;
}

void Geometry::findMassVector() {
    m_mass = Eigen::VectorXd::Zero(m_nn * m_nsd);

    // Total mass per small rectangle (2 triangular elements)
    m_mi = (m_rec_len * m_rec_wid * m_SimPar->thk() * m_SimPar->rho())
                / ((m_num_nodes_len - 1) * (m_num_nodes_wid - 1));

    for (int i = 0; i < m_num_nodes_len; i++) {
        for (int j = 0; j < m_num_nodes_wid; ++j) {

            // position in the mass vector
            int k = i + j * m_num_nodes_len;

            // corner nodes
            if ((i == 0 || i == m_num_nodes_len - 1) && (j == 0 || j == m_num_nodes_wid - 1)) {
                m_mass(k * m_nsd) = m_mi / 4.0;
                m_mass(k * m_nsd + 1) = m_mi / 4.0;
                m_mass(k * m_nsd + 2) = m_mi / 4.0;
            }
            // edge nodes
            else if ((i == 0 || i == m_num_nodes_len - 1) || (j == 0 || j == m_num_nodes_wid - 1)) {
                m_mass(k * m_nsd) = m_mi / 2.0;
                m_mass(k * m_nsd + 1) = m_mi / 2.0;
                m_mass(k * m_nsd + 2) = m_mi / 2.0;
            }
            // middle nodes
            else {
                m_mass(k * m_nsd) = m_mi;
                m_mass(k * m_nsd + 1) = m_mi;
                m_mass(k * m_nsd + 2) = m_mi;
            }
        }
    }
    /*
    // triangular element mass
    double tri_mass = (m_rec_len * m_rec_wid * m_SimPar->thk() * m_SimPar->rho())
           / (2.0 * (m_num_nodes_len - 1) * (m_num_nodes_wid - 1));
    m_mi = tri_mass / 3.0;

    for (int i = 0; i < m_nel; i++) {
        for (int j = 0; j < m_ndof; j++) {
            int inn = m_conn(i, j);

            m_mass((inn-1) * m_ndof) += m_mi;
            m_mass((inn-1) * m_ndof + 1) += m_mi;
            m_mass((inn-1) * m_ndof + 2) += m_mi;
        }
    }
    */

}

// print geometric information and generate connectivity file
void Geometry::printGeo() {
    /*
    std::cout << "----List of Nodes----" << '\n';
    for (int i = 0; i < m_nn; i++) {
        std::cout << i+1 << '\t';
        std::cout << std::setprecision(8) << std::fixed << m_dof(3*i) << '\t' << m_dof(3*i+1) << '\t' << m_dof(3*i+2) << std::endl;
    }
    std::cout << "----DOF Vector----" << '\n';
    std::cout << m_dof << std::endl;
    std::cout << "----Map of Nodes----" << '\n';
    std::cout << map_nodes << std::endl;
    std::cout << "----List of Elements----" << '\n';*/
    std::string filepath = m_SimPar->outputPath();
    std::string filename = "connectivity.txt";
    std::ofstream myfile((filepath+filename).c_str());

    for (int i = 0; i < m_nel; i++) {
        myfile << i+1 << '\t';
        for (int j = 0; j < m_nen; j++) {
            myfile << m_mesh[i][j] << '\t';
        }
        myfile << std::endl;
    }
}