#ifndef PLATES_SHELLS_ELEMENT_H
#define PLATES_SHELLS_ELEMENT_H

#include <vector>
#include <Eigen/Dense>
#include "node.h"
#include "parameters.h"

class Node;

class Element {
public:

    // constructor
    Element();
    Element(const unsigned int& num_el, Node* n1, Node* n2, Node* n3);

    // modifier
    void calculate_vec_edge();
    void calculate_len_edge();
    void calculate_area();
    void calculate_normal();
    void find_nearby_element(const Parameters& Params);

    // accessor
    unsigned int get_node_num(const int& num) const;
    double get_area() const;
    double get_len_edge(const int& num) const;
    Node* get_node(const int& num) const;
    Eigen::Vector3d get_vec_edge(const int& num) const;
    Eigen::Vector3d get_normal() const;

private:
    unsigned int m_num_el;
    unsigned int m_num_n1, m_num_n2, m_num_n3;
    int m_adj_element[3];
    double m_len_edge[3];
    double m_area;
    Eigen::Vector3d m_normal, m_dir12, m_dir23, m_dir13;

    Node* node1;
    Node* node2;
    Node* node3;
};

#endif //PLATES_SHELLS_ELEMENT_H
