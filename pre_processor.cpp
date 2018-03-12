#include "pre_processor.h"
#include "parameters.h"
#include "geometry.h"


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

void read_input(Parameters& sims) {
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
                sims.set_domain_len(std::stod(value_var));
            }
            else if (name_var == "rec_wid") {
                sims.set_domain_wid(std::stod(value_var));
            }
            else if (name_var == "num_nodes_len") {
                sims.set_num_nodes_len(std::stoul(value_var));
            }
            else if (name_var == "num_nodes_wid") {
                sims.set_num_nodes_wid(std::stoul(value_var));
            }
            else if (name_var == "dof") {
                sims.set_dof(std::stoi(value_var));
            }
            else if (name_var == "nen") {
                sims.set_nen(std::stoi(value_var));
            }
        }
    }

    std::cout << "Finish reading input file" << std::endl;
    sims.print_parameters();
}

// main pre-processor function
void pre_processor() {
    std::cout << "Pre-processor starts" << '\n';

    Parameters SimPar;
    read_input(SimPar);

    unsigned int num_nodes = SimPar.num_node_len() * SimPar.num_node_wid();
    unsigned int num_elements = 2 * (SimPar.num_node_len()-1) * (SimPar.num_node_wid()-1);
    int num_dof = SimPar.get_dof();
    int num_nen = SimPar.get_nen();

    Geometry SimGeo(num_nodes, num_elements, num_dof, num_nen);
    SimGeo.init_lst_nodes(SimPar);
    SimGeo.init_lst_elements(SimPar);
    SimGeo.print_geo(SimPar);
}