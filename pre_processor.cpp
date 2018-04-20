#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "pre_processor.h"
#include "parameters.h"
#include "global.h"


// ============================ //
//      Helper functions        //
// ============================ //

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

// ============================ //
//          Subroutines         //
// ============================ //

// read input file
void read_input(Parameters& Params) {

    std::string filepath = "/Users/chenjingyu/git/plates-shells/";
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

            if (name_var == "rec_len") {
                Params.set_rec_len(std::stod(value_var));             // length of the rectangular domain
            }
            else if (name_var == "rec_wid") {
                Params.set_rec_wid(std::stod(value_var));             // width of the rectangular domain
            }
            else if (name_var == "num_nodes_len") {
                Params.set_num_nodes_len(std::stoul(value_var));      // number of nodes along the length
            }
            else if (name_var == "num_nodes_wid") {
                Params.set_num_nodes_wid(std::stoul(value_var));      // number of nodes along the width
            }
            else if (name_var == "dof") {
                Params.set_ndof(std::stoi(value_var));                // degree of freedom
            }
            else if (name_var == "nen") {
                Params.set_nen(stoi(value_var));                      // number of nodes per element
            }
        }
    }

    Params.set_nn(Params.num_nodes_len() * Params.num_nodes_wid());
    Params.set_nel(2 * (Params.num_nodes_len()-1) * (Params.num_nodes_wid()-1));

    input_file.close();
}

// initialize coordinates (seeding)
void init_coord(Parameters& Params) {
    m_coord = Eigen::MatrixXd::Zero(Params.nn(), Params.ndof());
    m_dof = Eigen::VectorXd::Zero(Params.nn() * Params.ndof());
    map_nodes = Eigen::MatrixXi::Zero(Params.num_nodes_wid(), Params.num_nodes_len());

    double delta_len, delta_wid;
    delta_len = Params.rec_len() / (double) (Params.num_nodes_len() - 1);
    delta_wid = Params.rec_wid() / (double) (Params.num_nodes_wid() - 1);

    unsigned int node_count = 0;
    for (int i = 0; i < Params.num_nodes_wid(); i++) {
        for (int j = 0; j < Params.num_nodes_len(); j++) {
            if (node_count == Params.nn()) {
                break;
            }
            m_coord(node_count, 0) = (double) j * delta_len;
            m_coord(node_count, 1) = (double) i * delta_wid;
            m_coord(node_count, 2) = 0.0;

            m_dof(node_count * 3) = m_coord(node_count, 0);
            m_dof(node_count * 3 + 1) = m_coord(node_count, 1);
            m_dof(node_count * 3 + 2) = m_coord(node_count, 2);

            map_nodes(i,j) = node_count + 1;
            node_count++;
        }
    }
}

// initialize connectivity (meshing)
void init_conn(Parameters& Params) {
    m_conn = Eigen::MatrixXi::Zero(Params.nel(), Params.nen());

    unsigned int element_count = 0;
    int local_n1, local_n2, local_n3;
    int xpos, ypos;
    for (int i = 0; i < (Params.num_nodes_wid() - 1) * 2; i++) {
        for (int j = 0; j < (Params.num_nodes_len() - 1); j++) {
            if (i%2 == 0) {         // odd row
                ypos = j;
                xpos = (i+1)/2 + (i+1)%2 - 1;
                local_n1 = map_nodes(xpos,ypos);
                local_n2 = local_n1 + 1;
                local_n3 = local_n2 + Params.num_nodes_len();
            }
            else {                  // even row
                ypos = j + 1;
                xpos = (i+1)/2 + (i+1)%2;
                local_n1 = map_nodes(xpos,ypos);
                local_n2 = local_n1 - 1;
                local_n3 = local_n2 - Params.num_nodes_len();
            }
            m_conn(element_count, 0) = local_n1;
            m_conn(element_count, 1) = local_n2;
            m_conn(element_count, 2) = local_n3;
            element_count++;
        }
    }
}

// print geometric information
void print_geo(Parameters& Params) {
    std::cout << "----List of Nodes----" << '\n';
    for (int i = 0; i < Params.nn(); i++) {
        std::cout << i+1 << '\t';
        for (int j = 0; j < Params.ndof(); j++) {
            std::cout << std::setprecision(6) << std::fixed << m_coord(i,j) << '\t';
        }
        std::cout << std::endl;
    }
    std::cout << "----DOF Vector----" << '\n';
    std::cout << m_dof << std::endl;
    std::cout << "----Map of Nodes----" << '\n';
    std::cout << map_nodes << std::endl;
    std::cout << "----List of Elements----" << '\n';
    for (int i = 0; i < Params.nel(); i++) {
        std::cout << i+1 << '\t';
        for (int j = 0; j < Params.nen(); j++) {
            std::cout << m_conn(i,j) << '\t';
        }
        std::cout << std::endl;
    }
}

// ============================ //
//        Main functions        //
// ============================ //

// main pre-processor function
void pre_processor(Parameters& Params) {
    std::cout << "Pre-processor starts" << '\n';

    read_input(Params);
    std::cout << "Finish reading input file" << std::endl;
    Params.print_parameters();

    init_coord(Params);
    init_conn(Params);
    print_geo(Params);
}