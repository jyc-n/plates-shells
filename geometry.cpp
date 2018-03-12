#include "geometry.h"

Nodes::Nodes() {
    num_n_ = 0;
    x_ = 0.0;
    y_ = 0.0;
    z_ = 0.0;
}

void Nodes::set_num(const unsigned int &num) {
    num_n_ = num;
}

void Nodes::set_coord(const Eigen::Vector3d &coord) {
    x_ = coord(0);
    y_ = coord(1);
    z_ = coord(2);
}

Eigen::Vector3d& Nodes::get_coord() {
    Eigen::Vector3d coord;
    coord(0) = x_;
    coord(1) = y_;
    coord(2) = z_;
    return coord;
}

void Element::config_element(const unsigned int &num, std::vector<Nodes> lst) {
    node1.set_num(num_n1);
    node2.set_num(num_n2);
    node3.set_num(num_n3);

    node1.set_coord(lst[num_n1-1].get_coord());
    node2.set_coord(lst[num_n2-1].get_coord());
    node3.set_coord(lst[num_n3-1].get_coord());
}

 void Geometry::init_lst_nodes(Parameters& sims) {
     lst_nodes = Eigen::MatrixXd::Zero(nn_, dof_);

     double delta_len, delta_wid;
     delta_len = sims.rec_len() / (double) (sims.num_node_len() - 1);
     delta_wid = sims.rec_wid() / (double) (sims.num_node_wid() - 1);

     unsigned int node_count = 0;
     for (int i = 0; i < sims.num_node_wid(); i++) {
         for (int j = 0; j < sims.num_node_len(); j++) {
             if (node_count == nn_) {
                 break;
             }
             lst_nodes(node_count, 0) = (double) j * delta_len;
             lst_nodes(node_count, 1) = (double) i * delta_wid;
             lst_nodes(node_count, 2) = 0.0;
             node_count++;
         }
     }
 }

void Geometry::init_lst_elements(Parameters &sims) {
    //
}

void Geometry::print_nodes() {
    std::cout << "----List of Nodes----" << '\n';
    for (int i = 0; i < nn_; i++) {
        std::cout << i << '\t';
        for (int j = 0; j < dof_; j++) {
            std::cout << std::setprecision(6) << std::fixed << lst_nodes(i,j) << '\t';
        }
        std::cout << std::endl;
    }
}

