#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cassert>
#include "pre_processor.h"
#include "parameters.h"
#include "geometry.h"
#include "node.h"
#include "element.h"
#include "hinge.h"
#include "edge.h"

// ========================================= //
//      Declaration of helper functions      //
// ========================================= //

std::string& trim(std::string& str);


// ========================================= //
//             Member functions              //
// ========================================= //

PreProcessorImpl::PreProcessorImpl(Parameters* SimPar, Geometry* SimGeo) {
    m_SimPar = SimPar;
    m_SimGeo = SimGeo;
}

// main pre-processor function
void PreProcessorImpl::PreProcess() {
    //std::cout << "Pre-processor starts" << '\n';

    readInput();
    //std::cout << "Finish reading input file" << std::endl;
    m_SimPar->print_parameters();

    m_SimGeo->buildGeo();
    initCoord();
    initConn();
    initGeoList();
    m_SimGeo->calcMass();
    m_SimGeo->printGeo();
}

// read input file
void PreProcessorImpl::readInput() {

    std::string filepath = m_SimPar->inputPath();
    std::string fileid = "input.txt";

    std::ifstream input_file((filepath+fileid).c_str());

    if (!input_file.good()) {
        std::cout << "Error opening " << fileid << '\n';
        throw "Input file name error";
    }

    while (!input_file.eof()) {
        std::string line;
        getline(input_file,line);
        if (line.length() != 0) {
            trim(line);

            std::string name_var = line.substr(0,line.find('='));
            std::string value_var = line.substr(line.find('=')+1);

            if (name_var == "rec_len")
                m_SimGeo->set_rec_len(std::stod(value_var));             // length of the rectangular domain
            else if (name_var == "rec_wid")
                m_SimGeo->set_rec_wid(std::stod(value_var));             // width of the rectangular domain
            else if (name_var == "num_nodes_len")
                m_SimGeo->set_num_nodes_len(std::stoul(value_var));      // number of nodes along the length
            else if (name_var == "num_nodes_wid")
                m_SimGeo->set_num_nodes_wid(std::stoul(value_var));      // number of nodes along the width
            else if (name_var == "dof")
                m_SimGeo->set_ndof(std::stoi(value_var));                // degree of freedom
            else if (name_var == "nen")
                m_SimGeo->set_nen(std::stoi(value_var));                 // number of nodes per element
            else if (name_var == "dt")
                m_SimPar->set_dt(std::stod(value_var));                  // step size
            else if (name_var == "nst")
                m_SimPar->set_nst(std::stoi(value_var));                 // total number of steps
            else if (name_var == "iter_lim")
                m_SimPar->set_iter_lim(std::stoi(value_var));            // maximum number of iterations allowed per time step
            else if (name_var == "ctol")
                m_SimPar->set_ctol(std::stod(value_var));                // tolerance scaling factor
            else if (name_var == "E_modulus")
                m_SimPar->set_E_modulus(std::stod(value_var));           // Young's modulus
            else if (name_var == "nu")
                m_SimPar->set_nu(std::stod(value_var));                  // Poisson's ratio
            else if (name_var == "rho")
                m_SimPar->set_rho(std::stod(value_var));                 // density
            else if (name_var == "thk")
                m_SimPar->set_thk(std::stod(value_var));                 // thickness
            else if (name_var == "vis")
                m_SimPar->set_vis(std::stod(value_var));                 // viscosity
            else if (name_var == "gconst")
                m_SimPar->set_gconst(std::stod(value_var));              // gravitational acceleration constant
            else if (name_var == "outop")
                m_SimPar->set_outop((bool) std::stoi(value_var));        // output option
            else if (name_var == "solver_op")
                m_SimPar->set_solver_op((bool) std::stoi(value_var));    // solver option
        }
    }
    m_SimPar->set_kstretch();
    m_SimPar->set_kshear();
    m_SimPar->set_kbend();

    m_SimGeo->set_nn(m_SimGeo->num_nodes_len() * m_SimGeo->num_nodes_wid());
    m_SimGeo->set_nel(2 * (m_SimGeo->num_nodes_len()-1) * (m_SimGeo->num_nodes_wid()-1));
    m_SimGeo->set_nhinge((3 * m_SimGeo->nel() - 2 * (m_SimGeo->num_nodes_len() + m_SimGeo->num_nodes_wid() - 2)) / 2);
    m_SimGeo->set_nedge(3 * m_SimGeo->nel() - m_SimGeo->nhinge());

    input_file.close();
}

// initialize coordinates (seeding)
void PreProcessorImpl::initCoord() {

    double delta_len, delta_wid;
    delta_len = m_SimGeo->rec_len() / (double) (m_SimGeo->num_nodes_len() - 1);
    delta_wid = m_SimGeo->rec_wid() / (double) (m_SimGeo->num_nodes_wid() - 1);

    unsigned int node_count = 0;
    for (int i = 0; i < m_SimGeo->num_nodes_wid(); i++) {
        for (int j = 0; j < m_SimGeo->num_nodes_len(); j++) {

            if (node_count == m_SimGeo->nn())
                break;

            m_SimGeo->m_coord(node_count, 0) = (double) j * delta_len;
            m_SimGeo->m_coord(node_count, 1) = (double) i * delta_wid;
            m_SimGeo->m_coord(node_count, 2) = 0.0;

            m_SimGeo->m_dof(node_count * 3)     = m_SimGeo->m_coord(node_count, 0);
            m_SimGeo->m_dof(node_count * 3 + 1) = m_SimGeo->m_coord(node_count, 1);
            m_SimGeo->m_dof(node_count * 3 + 2) = m_SimGeo->m_coord(node_count, 2);

            m_SimGeo->map_nodes(i,j) = node_count + 1;
            node_count++;
        }
    }
}

// initialize connectivity (meshing)
void PreProcessorImpl::initConn() {

    unsigned int element_count = 0;
    int local_n1, local_n2, local_n3;
    int xpos, ypos;
    for (int i = 0; i < (m_SimGeo->num_nodes_wid() - 1) * 2; i++) {
        for (int j = 0; j < (m_SimGeo->num_nodes_len() - 1); j++) {
            if (i%2 == 0) {         // odd row
                ypos = j;
                xpos = (i+1)/2 + (i+1)%2 - 1;
                local_n1 = m_SimGeo->map_nodes(xpos,ypos);
                local_n2 = local_n1 + 1;
                local_n3 = local_n2 + m_SimGeo->num_nodes_len();
            }
            else {                  // even row
                ypos = j + 1;
                xpos = (i+1)/2 + (i+1)%2;
                local_n1 = m_SimGeo->map_nodes(xpos,ypos);
                local_n2 = local_n1 - 1;
                local_n3 = local_n2 - m_SimGeo->num_nodes_len();
            }
            m_SimGeo->m_conn(element_count, 0) = local_n1;
            m_SimGeo->m_conn(element_count, 1) = local_n2;
            m_SimGeo->m_conn(element_count, 2) = local_n3;
            element_count++;
        }
    }
}

// initialize element list
void PreProcessorImpl::initGeoList() {

    // build node list
    for (unsigned int i = 0; i < m_SimGeo->nn(); i++) {
        Node temp(i+1, m_SimGeo->m_coord(i,0), m_SimGeo->m_coord(i,1), m_SimGeo->m_coord(i,2));
        m_SimGeo->m_nodeList.push_back(temp);
    }
    assert(m_SimGeo->m_nodeList.size() == m_SimGeo->nn());

    // build element list
    for (unsigned int i = 0; i < m_SimGeo->nel(); i++) {
        // pointers to 3 node object
        Node* pn1 = &m_SimGeo->m_nodeList[m_SimGeo->m_conn(i,0)-1];
        Node* pn2 = &m_SimGeo->m_nodeList[m_SimGeo->m_conn(i,1)-1];
        Node* pn3 = &m_SimGeo->m_nodeList[m_SimGeo->m_conn(i,2)-1];

        Element temp(i+1, pn1, pn2, pn3);

        temp.find_nearby_element(*m_SimGeo);
        m_SimGeo->m_elementList.push_back(temp);
    }
    assert(m_SimGeo->m_elementList.size() == m_SimGeo->nel());

    // find all edges and hinges
    for (unsigned int i = 0; i < m_SimGeo->nel(); i++) {
        m_SimGeo->m_elementList[i].find_edges(m_SimGeo->m_elementList);
        m_SimGeo->m_elementList[i].find_hinges(m_SimGeo->m_elementList);
    }

    // build edge list
    for (unsigned int i = 0; i < m_SimGeo->nel(); i++) {
        for (int j = 1; j <= 3; j++) {
            Edge* tempEdge = m_SimGeo->m_elementList[i].get_edge(j);
            if (tempEdge != nullptr && !tempEdge->check_visited()) {
                tempEdge->mark_visited();
                m_SimGeo->m_edgeList.push_back(tempEdge);
            }

            Hinge* tempHinge = m_SimGeo->m_elementList[i].get_hinge(j);
            if (tempHinge != nullptr && !tempHinge->check_visited()) {
                tempHinge->mark_visited();
                m_SimGeo->m_hingeList.push_back(tempHinge);
            }
        }
    }
    assert(m_SimGeo->edgeNumCheck());
    assert(m_SimGeo->hingeNumCheck());
}

// ========================================= //
//     Implementation of helper functions    //
// ========================================= //

// clear all unnecessary characters in the string
std::string& trim(std::string& str) {
    if (str.empty()) {
        return str;
    }
    str.erase(str.begin() + str.find('!'), str.end());
    str.erase(str.find_last_not_of(' ')+1);
    auto pos = str.find(' ');
    while (pos != std::string::npos) {
        str.replace(pos, 1, "");
        pos = str.find(' ');
    }
    return str;
}