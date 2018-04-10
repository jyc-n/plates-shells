#include <iostream>
#include <iomanip>
#include "geometry.h"

void Geometry::init_lst_nodes(Parameters& Params) {
    lst_nodes = Eigen::MatrixXd::Zero(Params.nn(), Params.ndof());
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
            lst_nodes(node_count, 0) = (double) j * delta_len;
            lst_nodes(node_count, 1) = (double) i * delta_wid;
            lst_nodes(node_count, 2) = 0.0;
            map_nodes(i,j) = node_count + 1;
            node_count++;
        }
    }
}

void Geometry::init_lst_elements(Parameters& Params) {
    lst_elements = Eigen::MatrixXi::Zero(Params.nel(), Params.nen());

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
            lst_elements(element_count, 0) = local_n1;
            lst_elements(element_count, 1) = local_n2;
            lst_elements(element_count, 2) = local_n3;
            element_count++;
        }
    }
}

Eigen::MatrixXd Geometry::lst_coord() const {
    return lst_nodes;
}

Eigen::MatrixXi Geometry::lst_conn() const {
    return lst_elements;
}

void Geometry::print_geo(Parameters& Params) {
    std::cout << "----List of Nodes----" << '\n';
    for (int i = 0; i < Params.nn(); i++) {
        std::cout << i+1 << '\t';
        for (int j = 0; j < Params.ndof(); j++) {
            std::cout << std::setprecision(6) << std::fixed << lst_nodes(i,j) << '\t';
        }
        std::cout << std::endl;
    }
    std::cout << "----Map of Nodes----" << '\n';
    std::cout << map_nodes << std::endl;
    std::cout << "----List of Elements----" << '\n';
    for (int i = 0; i < Params.nel(); i++) {
        std::cout << i+1 << '\t';
        for (int j = 0; j < Params.nen(); j++) {
            std::cout << lst_elements(i,j) << '\t';
        }
        std::cout << std::endl;
    }
}