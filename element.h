#ifndef PLATES_SHELLS_ELEMENT_H
#define PLATES_SHELLS_ELEMENT_H

#include <Eigen/Dense>
#include "geometry.h"

class Node;
class Hinge;

class Element {
public:

    // constructor
    Element();
    Element(const unsigned int num_el, Node* n1, Node* n2, Node* n3);

    // modifier
    void calculate_vec_edge();
    void calculate_len_edge();
    void calculate_area();
    void calculate_normal();
    void find_nearby_element(Geometry& Geo);
    void find_hinges(Element* l_element);

    // accessor
    unsigned int get_node_num(const int num) const;
    double get_area() const;
    double get_len_edge(const int num) const;
    Node* get_node(const int num) const;
    int get_nearby_element(int num) const;
    int get_overlapped_edge_num(int num_el) const;
    int get_remain_node_num(int num_edge) const;
    bool is_hinge(int num_edge) const;
    Hinge* get_hinge(int num_edge) const;
    Eigen::Vector3d get_vec_edge(const int num) const;
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
    Hinge* m_hinges[3];
};

#endif //PLATES_SHELLS_ELEMENT_H
