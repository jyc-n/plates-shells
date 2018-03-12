#ifndef PLATES_SHELLS_GEOMETRY_H
#define PLATES_SHELLS_GEOMETRY_H

#include <iomanip>
#include <vector>
#include <Eigen/Dense>

#include "parameters.h"

class Nodes {
public:
    Nodes();
    /*
    Nodes(const unsigned int& num_n,
          const double& xpos,
          const double& ypos,
          const double& zpos)
            : num_n_(num_n), x_(xpos), y_(ypos), z_(zpos) {}
            */
    void set_num(const unsigned int& num);
    void set_coord(const Eigen::Vector3d& coord);
    Eigen::Vector3d& get_coord();

private:
    unsigned int num_n_;
    double x_, y_, z_;
};


class Element {
public:
    Element(const unsigned int& num_el,
            const unsigned int& nn1,
            const unsigned int& nn2,
            const unsigned int& nn3)
             : num_el_(num_el), num_n1(nn1), num_n2(nn2), num_n3(nn3) {}
    void config_element(const unsigned int& num, std::vector<Nodes> lst);

private:
    unsigned int num_el_;
    unsigned int num_n1, num_n2, num_n3;
    Nodes node1, node2, node3;
};


class Geometry {
public:
    Geometry(const unsigned int& nn, const unsigned int& nel, const int& dof, const int& nen)
            : nn_(nn), nel_(nel), ndof_(dof), nen_(nen) {}

    void init_lst_nodes(Parameters& sims);
    void init_lst_elements(Parameters& sims);
    void print_geo(Parameters& sims);

private:
    unsigned int nn_;    // number of nodes
    unsigned int nel_;   // number of elements
    int ndof_;           // degree of freedom
    int nen_;            // number of nodes per element

    Eigen::MatrixXi lst_elements;
    Eigen::MatrixXd lst_nodes;
    Eigen::MatrixXi map_nodes;
};

#endif //PLATES_SHELLS_GEOMETRY_H
