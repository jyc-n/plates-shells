#include "parameters.h"

Parameters::Parameters() {
    rec_len_ = 0.0;
    rec_wid_ = 0.0;
    num_nodes_len_ = 0;
    num_nodes_wid_ = 0;
}

void Parameters::print_parameters() {
    std::cout << "List of parameters" << '\n';
    std::cout << "length of the rectangular domain: " << rec_len() << '\n';
    std::cout << "width of the rectangular domain: " << rec_wid() << '\n';
    std::cout << "number of nodes along the length: " << num_node_len() << '\n';
    std::cout << "number of nodes along the width: " << num_node_wid() << '\n';
}