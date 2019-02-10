#ifndef PLATES_SHELLS_ELEMENT_H
#define PLATES_SHELLS_ELEMENT_H

#include <vector>
#include <Eigen/Dense>
#include "geometry.h"

class Node;
class Edge;
class Hinge;

class Element {
public:

    // constructor
    Element(const unsigned int num_el, Node* n1, Node* n2, Node* n3);
    ~Element();

    // modifier
    void calculate_phi0();
    void calculate_area();
    void find_edge_index();
    void find_adjacent_element(const int index, const int num_element);
    Edge* build_edges(const int index);
    Hinge* build_hinges(const int index, Element* adj_element);

    // accessor
    unsigned int get_node_num(int num) const;
    int get_edge_index(int num) const;
    int get_element_num() const;
    double get_phi0() const;
    double get_area() const;
    Node* get_node(int num) const;

    double m_k;

private:
    int get_which_edge(int num) const;
    int get_remain_node_num(int num_edge) const;

    unsigned int m_num_el;
    int m_adj_element[3];
    int m_edgeIndex[3];
    double m_phi0;
    double m_area;

    Node* m_node1;
    Node* m_node2;
    Node* m_node3;
};

inline int Element::get_element_num() const { return m_num_el; }
inline double Element::get_phi0() const { return m_phi0; }
inline double Element::get_area() const { return m_area; }

#endif //PLATES_SHELLS_ELEMENT_H
