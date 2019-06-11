#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <string>
#include <cassert>
#include <unordered_map>
#include "pre_processor.h"
#include "arguments.h"
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
void PreProcessorImpl::PreProcess(Arguments t_args) {
    // read input file
    readInput();
    // specify number of nodes
    if (t_args.argc == NUMS) {
        // m_SimGeo->set_num_nodes_len(t_args.num_len+1);      //! only used when edge is clamped
        m_SimGeo->set_num_nodes_len(t_args.num_len);
        m_SimGeo->set_num_nodes_wid(t_args.num_wid);
    }
    // specify number of nodes and dimensions
    else if (t_args.argc == NUMS_DIMS) {
        m_SimGeo->set_num_nodes_len(t_args.num_len);
        m_SimGeo->set_num_nodes_wid(t_args.num_wid);
        m_SimGeo->set_rec_len(t_args.len);
        m_SimGeo->set_rec_wid(t_args.wid);
    }
    // m_SimPar->find_fullOutputPath(m_SimGeo->num_nodes_len()-1, m_SimGeo->num_nodes_wid()); // ! only -1 when clamped
    m_SimPar->find_fullOutputPath(m_SimGeo->num_nodes_len(), m_SimGeo->num_nodes_wid());

    m_SimGeo->set_nn();
    m_SimGeo->set_nel();
    m_SimGeo->set_nhinge();
    m_SimGeo->set_nedge();

    print_parameters();

    std::cout << "input file read, preprocessor starts" << std::endl;
    buildNodes(t_args.argc);
    std::cout << "seeding completed" << std::endl;
    buildMesh();
    std::cout << "meshing completed\nbuilding node, element, edge, hinge lists (this may take a while)" << std::endl;
    buildNodeElementList();
    buildEdgeHingeList();
    std::cout << "all lists building completed" << std::endl;
    m_SimGeo->findMassVector();
    std::cout << "mass calculated completed" << std::endl;
    m_SimGeo->writeConnectivity();
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
            if (line.empty())
                continue;

            std::string name_var = line.substr(0,line.find('='));
            std::string value_var = line.substr(line.find('=')+1);

            // TODO: use switch!! hash the string

            if (name_var == "rec_len")
                m_SimGeo->set_rec_len(std::stod(value_var));             // length of the rectangular domain
            else if (name_var == "rec_wid")
                m_SimGeo->set_rec_wid(std::stod(value_var));             // width of the rectangular domain
            else if (name_var == "num_nodes_len")
                m_SimGeo->set_num_nodes_len(std::stoul(value_var));      // number of nodes along the length
            else if (name_var == "num_nodes_wid")
                m_SimGeo->set_num_nodes_wid(std::stoul(value_var));      // number of nodes along the width
            else if (name_var == "dx")
                m_SimGeo->set_dx(std::stod(value_var));                  // spatial step size
            else if (name_var == "datum_plane")
                m_SimGeo->set_datum(std::stoi(value_var));               // datum plane
            else if (name_var == "angle")
                m_SimGeo->set_angle(std::stod(value_var));               // rotation about origin
            else if (name_var == "dof")
                m_SimGeo->set_nsd(std::stoi(value_var));                 // 2D or 3D
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

    input_file.close();
}

// read Abaqus .inp file
void PreProcessorImpl::readGeoFile() {
    // TODO: implement read abaqus input file subroutine
}

// initialize coordinates (seeding)
void PreProcessorImpl::buildNodes(int opt) {
    double dx1 = 0, dx2 = 0;

    // if (AR_FLAG) {
    //     dx1 = m_SimGeo->dx();
    //     dx2 = m_SimGeo->dx();
    // }
    // else {
    //     // increment along length/width
    //     dx1 = m_SimGeo->rec_len() / (double) (m_SimGeo->num_nodes_len() - 1);
    //     dx2 = m_SimGeo->rec_wid() / (double) (m_SimGeo->num_nodes_wid() - 1);
    // }
    // increment along length/width
    dx1 = m_SimGeo->rec_len() / (double) (m_SimGeo->num_nodes_len() - 1);
    dx2 = m_SimGeo->rec_wid() / (double) (m_SimGeo->num_nodes_wid() - 1);
    // ! only use the following for clamped case
    // dx1 = dx2;

    // position of datum plane
    int datum_op = m_SimGeo->datum();

    // rotation of the plane (2D rotation matrix)
    //    [ T1   -T2 ]  [ q1 ]   =   [ T1 * q1 - T2 * q2 ]
    //    [ T2    T1 ]  [ q2 ]       [ T2 * q1 + T1 * q2 ]
    double T1 = cos(m_SimGeo->angle() * M_PI / 180.0);
    double T2 = sin(m_SimGeo->angle() * M_PI / 180.0);

    m_SimGeo->m_nodes.resize(m_SimGeo->nn());
    for (int i = 0; i < m_SimGeo->nn(); i++)
        m_SimGeo->m_nodes[i].fill(0.0);

    for (int i = 0; i < m_SimGeo->num_nodes_wid(); i++) {
        for (int j = 0; j < m_SimGeo->num_nodes_len(); j++) {
            // coordinates in rotated plane
            double q1 = (double) j * dx1;
            double q2 = (double) i * dx2;

            // coordinates in physical plane
            double x = 0.0, y = 0.0, z = 0.0;
            switch (datum_op) {
                default:
                case 1:             // xy plane
                    x = T1 * q1 - T2 * q2;
                    y = T2 * q1 + T1 * q2;
                    z = 0.0;
                    break;
                case 2:             // xz plane
                    x = T1 * q1 - T2 * q2;
                    y = 0.0;
                    z = T2 * q1 + T1 * q2;
                    break;
                case 3:             // yz plane
                    x = 0.0;
                    y = T1 * q1 - T2 * q2;
                    z = T2 * q1 + T1 * q2;
                    break;
            }
            int index = i * m_SimGeo->num_nodes_len() + j;
            assert(index >= 0 && index <= m_SimGeo->nn());
            m_SimGeo->m_nodes[index][0] = x;
            m_SimGeo->m_nodes[index][1] = y;
            m_SimGeo->m_nodes[index][2] = z;
        }
    }
    // m_SimGeo->translateNodes(0, -dx1); //! only used when clamped!
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

// initialize node and element lists
void PreProcessorImpl::buildNodeElementList() {

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
        m_SimGeo->m_elementList.emplace_back(temp);
    }
    assert(m_SimGeo->m_elementList.size() == m_SimGeo->nel());
}

// initialize edge and hinge lists
void PreProcessorImpl::buildEdgeHingeList() {
    std::unordered_multimap<int, int> edgeHashTable;
    std::vector<int> edgeIndex;
    // use open hash table to store edge indices and the adjacent elements
    // use an additional vector to store all the keys (edge indices) of the hash table
    for (auto &iel : m_SimGeo->m_elementList) {
        // find and save all edge indices
        iel.find_edge_index();
        for (int iedge = 0; iedge < 3; iedge++) {
            edgeHashTable.emplace(iel.get_edge_index(iedge), iel.get_element_num());
            edgeIndex.push_back(iel.get_edge_index(iedge));
        }
        /*
        std::cout << (*iel).get_element_num() << '\n';
        std::cout << (*iel).get_node_num(1) << '\t' << (*iel).get_node_num(2) << '\t' << (*iel).get_edge_index(0) << '\n';
        std::cout << (*iel).get_node_num(2) << '\t' << (*iel).get_node_num(3) << '\t' << (*iel).get_edge_index(1) << '\n';
        std::cout << (*iel).get_node_num(1) << '\t' << (*iel).get_node_num(3) << '\t' << (*iel).get_edge_index(2) << "\n\n";
        */
    }
    // remove the repeated keys
    std::sort(edgeIndex.begin(), edgeIndex.end());
    auto last = std::unique(edgeIndex.begin(), edgeIndex.end());
    edgeIndex.erase(last, edgeIndex.end());

    // traverse the non-empty buckets of the hash table
    for (int &key : edgeIndex) {
        //std::cout << "key is " << *key << '\n';

        // store {#element1, #element2}
        Eigen::Vector2i edgeInfo;
        edgeInfo.fill(-1);

        // locate the bucket with given key
        auto bucket_range = edgeHashTable.equal_range(key);
        int i = 0;
        for (auto p = bucket_range.first; p != bucket_range.second; p++, i++) {
            edgeInfo[i] = p->second;
            //std::cout << "Key:[" << p->first << "] Value:[" << p->second << "]\n";
        }

        // #element1 must have actual element number
        // #element2 can be -1 (which means the edge isn't shared by 2 elements. It's an edge)
        assert(edgeInfo[0] != -1);

        // first, build the edge and store the pointer to the list
        Element* this_element = &m_SimGeo->m_elementList[edgeInfo[0]-1];
        this_element->find_adjacent_element(key, edgeInfo[1]);
        m_SimGeo->m_edgeList.emplace_back(this_element->build_edges(key));

        // if this is a hinge, build the hinge and store the pointer to the list
        if (edgeInfo[1] != -1) {
            Element* adj_element = &m_SimGeo->m_elementList[edgeInfo[1]-1];
            adj_element->find_adjacent_element(key, edgeInfo[0]);
            m_SimGeo->m_hingeList.emplace_back(this_element->build_hinges(key, adj_element));
        }
    }
    // NOTE: edge/hinge number check only works for structured rectangular mesh
    assert(m_SimGeo->edgeNumCheck());
    assert(m_SimGeo->hingeNumCheck());
}

void PreProcessorImpl::print_parameters() {
    std::cout << "\n---------------------------------------" << std::endl;
    std::cout << "\t\tList of parameters" << '\n';
    std::cout << "\t\t# of Nodes:       " << m_SimGeo->nn() << '\n';
    std::cout << "\t\t# of Elements:    " << m_SimGeo->nel() << '\n';
    std::cout << "\t\tStep size:        " << m_SimPar->dt() << '\n';
    std::cout << "\t\tTotal # of steps: " << m_SimPar->nst() << '\n';
    std::cout << "---------------------------------------\n" << std::endl;
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