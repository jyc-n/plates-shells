#include <iostream>
#include <fstream>
#include <string>
#include "pre_processor.h"
#include "parameters.h"
#include "geometry.h"
#include "element.h"
#include "node.h"

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

// main pre-processor function
std::vector<Element> pre_processor(Geometry& Geo, Parameters& Params) {
    std::cout << "Pre-processor starts" << '\n';

    read_input(Params);
    std::cout << "Finish reading input file" << std::endl;
    Params.print_parameters();

    Geo.init_lst_nodes(Params);
    Geo.init_lst_elements(Params);
    Geo.print_geo(Params);

    std::vector<Node> init_node_list;
    for (int i = 0; i < Params.nn(); i++) {
        Node new_node(Geo.lst_coord()(i,0), Geo.lst_coord()(i,1), Geo.lst_coord()(i,2));
        init_node_list.push_back(new_node);
    }
    std::vector<Element> init_element_list;
    for (int i = 0; i < Params.nel(); i++) {
        Element new_element(i, Geo.lst_conn()(i,0), Geo.lst_conn()(i,1), Geo.lst_conn()(i,2));
        new_element.set_node(init_node_list);

        new_element.calculate_dir();
        new_element.calculate_angle();
        new_element.calculate_normal();

        init_element_list.push_back(new_element);
    }
    return init_element_list;
}