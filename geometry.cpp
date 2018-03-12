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
     lst_nodes = Eigen::MatrixXd::Zero(nn_, ndof_);
     map_nodes = Eigen::MatrixXi::Zero(sims.num_node_wid(), sims.num_node_len());

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
             map_nodes(i,j) = node_count + 1;
             node_count++;
         }
     }
 }

void Geometry::init_lst_elements(Parameters& sims) {
    lst_elements = Eigen::MatrixXi::Zero(nel_, nen_);

    unsigned int element_count = 0;
    int local_n1, local_n2, local_n3;
    int xpos, ypos;
    for (int i = 0; i < (sims.num_node_wid()-1)*2; i++) {
        for (int j = 0; j < (sims.num_node_len()-1); j++) {
            if (i%2 == 0) {         // odd row
                ypos = j;
                xpos = (i+1)/2 + (i+1)%2 - 1;
                local_n1 = map_nodes(xpos,ypos);
                local_n2 = local_n1 + 1;
                local_n3 = local_n2 + sims.num_node_len();
            }
            else {                  // even row
                ypos = j + 1;
                xpos = (i+1)/2 + (i+1)%2;
                local_n1 = map_nodes(xpos,ypos);
                local_n2 = local_n1 - 1;
                local_n3 = local_n2 - sims.num_node_len();
            }
            lst_elements(element_count, 0) = local_n1;
            lst_elements(element_count, 1) = local_n2;
            lst_elements(element_count, 2) = local_n3;
            element_count++;
        }
    }

}

void Geometry::print_geo(Parameters& sims) {
    std::cout << "----List of Nodes----" << '\n';
    for (int i = 0; i < nn_; i++) {
        std::cout << i+1 << '\t';
        for (int j = 0; j < ndof_; j++) {
            std::cout << std::setprecision(6) << std::fixed << lst_nodes(i,j) << '\t';
        }
        std::cout << std::endl;
    }
    std::cout << "----Map of Nodes----" << '\n';
    std::cout << map_nodes << std::endl;
    std::cout << "----List of Elements----" << '\n';
    for (int i = 0; i < nel_; i++) {
        std::cout << i+1 << '\t';
        for (int j = 0; j < nen_; j++) {
            std::cout << lst_elements(i,j) << '\t';
        }
        std::cout << std::endl;
    }
}

