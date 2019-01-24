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
    // read input file
    readInput();
    m_SimPar->print_parameters();

    std::cout << "input file read, preprocessor starts" << std::endl;
    buildNodes();
    std::cout << "seeding completed" << std::endl;
    buildMesh();
    std::cout << "meshing completed\nbuilding node, element, edge, hinge lists (this may take a while)" << std::endl;
    buildGeoList();
    std::cout << "all lists building completed" << std::endl;
    m_SimGeo->findMassVector();
    std::cout << "mass calculated completed" << std::endl;
    m_SimGeo->printGeo();
    std::cout << "connectivity file written, preprocessing completed\n" << std::endl;
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
            else if (name_var == "datum_plane")
                m_SimGeo->set_datum(std::stoi(value_var));               // degree of freedom
            else if (name_var == "dof")
                m_SimGeo->set_nsd(std::stoi(value_var));                // degree of freedom
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
            else if (name_var == "out_freq")
                m_SimPar->set_out_freq(std::stoi(value_var));            // output frequency
            else if (name_var == "solver_op")
                m_SimPar->set_solver_op((bool) std::stoi(value_var));    // solver option
            else if (name_var == "info_style")
                m_SimPar->set_info_style((bool) std::stoi(value_var));   // console output style
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

// read Abaqus .inp file
void PreProcessorImpl::readGeoFile() {
    // TODO: implement read abaqus input file subroutine
}

// initialize coordinates (seeding)
void PreProcessorImpl::buildNodes() {
    // increment along length/width
    double delta_len = m_SimGeo->rec_len() / (double) (m_SimGeo->num_nodes_len() - 1);
    double delta_wid = m_SimGeo->rec_wid() / (double) (m_SimGeo->num_nodes_wid() - 1);
    // position of datum plane
    int datum_op = m_SimGeo->datum();

    m_SimGeo->m_nodes.resize(m_SimGeo->nn());
    for (int i = 0; i < m_SimGeo->nn(); i++)
        m_SimGeo->m_nodes[i].fill(0.0);

    for (int i = 0; i < m_SimGeo->num_nodes_wid(); i++) {
        for (int j = 0; j < m_SimGeo->num_nodes_len(); j++) {
            double x = 0.0, y = 0.0, z = 0.0;
            switch (datum_op) {
                default:
                case 1:             // xy plane
                    x = (double) j * delta_len;
                    y = (double) i * delta_wid;
                    z = 0.0;
                    break;
                case 2:             // xz plane
                    x = (double) j * delta_len;
                    y = 0.0;
                    z = -(double) i * delta_wid;
                    break;
                case 3:             // yz plane
                    x = 0.0;
                    y = (double) j * delta_len;
                    z = -(double) i * delta_wid;
                    break;
            }
            int index = i * m_SimGeo->num_nodes_len() + j;
            assert(index >= 0 && index <= m_SimGeo->nn());
            m_SimGeo->m_nodes[index][0] = x;
            m_SimGeo->m_nodes[index][1] = y;
            m_SimGeo->m_nodes[index][2] = z;
        }
    }
    assert(m_SimGeo->m_nodes.size() == m_SimGeo->nn());
}

// initialize connectivity (meshing)
void PreProcessorImpl::buildMesh() {
    int num_wid = m_SimGeo->num_nodes_wid();
    int num_len = m_SimGeo->num_nodes_len();

    m_SimGeo->m_mesh.resize(m_SimGeo->nel());
    for (int i = 0; i < m_SimGeo->nel(); i++)
        m_SimGeo->m_mesh[i].fill((int)0);

    for (int i = 0; i < num_wid-1; i++) {
        for (int j = 0; j < num_len-1; j++) {
            int el1 = 2 * ((num_len-1) * i + j);
            int el2 = el1 + 1;
            assert(el1 < m_SimGeo->nel() && el2 < m_SimGeo->nel());
            m_SimGeo->m_mesh[el1] << i*num_len+j+1, i*num_len+j+2, (i+1)*num_len+j+2;
            m_SimGeo->m_mesh[el2] << (i+1)*num_len+j+2, (i+1)*num_len+j+1, i*num_len+j+1;
        }
    }
    assert(m_SimGeo->m_mesh.size() == m_SimGeo->nel());
}

// initialize element list
void PreProcessorImpl::buildGeoList() {

    // build node list
    for (unsigned int i = 0; i < m_SimGeo->nn(); i++) {
        Node temp(i+1, &m_SimGeo->m_nodes[i]);
        m_SimGeo->m_nodeList.emplace_back(temp);
    }
    assert(m_SimGeo->m_nodeList.size() == m_SimGeo->nn());

    // build element list
    for (unsigned int i = 0; i < m_SimGeo->nel(); i++) {
        // pointers to 3 node object
        Node* pn1 = &m_SimGeo->m_nodeList[m_SimGeo->m_mesh[i][0]-1];
        Node* pn2 = &m_SimGeo->m_nodeList[m_SimGeo->m_mesh[i][1]-1];
        Node* pn3 = &m_SimGeo->m_nodeList[m_SimGeo->m_mesh[i][2]-1];

        Element temp(i+1, pn1, pn2, pn3);

        temp.find_nearby_element(m_SimGeo->nel(), m_SimGeo->m_mesh);
        m_SimGeo->m_elementList.emplace_back(temp);
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
                m_SimGeo->m_edgeList.emplace_back(tempEdge);
            }

            Hinge* tempHinge = m_SimGeo->m_elementList[i].get_hinge(j);
            if (tempHinge != nullptr && !tempHinge->check_visited()) {
                tempHinge->mark_visited();
                m_SimGeo->m_hingeList.emplace_back(tempHinge);
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