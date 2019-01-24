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
    Element();
    Element(const unsigned int num_el, Node* n1, Node* n2, Node* n3);
    ~Element();

    // modifier
    void calculate_phi0();
    void calculate_area();
    void find_nearby_element(int nel, const VectorMesh& mesh);
    void find_edges(std::vector<Element>& l_element);
    void find_hinges(std::vector<Element>& l_element);

    // accessor
    unsigned int get_node_num(const int num) const;
    int get_nearby_element(int num) const;
    bool is_hinge(int num_edge) const;
    double get_phi0() const;
    double get_area() const;
    Node* get_node(const int num) const;
    Edge* get_edge(int num_edge) const;
    Hinge* get_hinge(int num_edge) const;

    double m_k;

private:

    int get_overlapped_edge_num(int num_el) const;
    int get_remain_node_num(int num_edge) const;

    unsigned int m_num_el;
    unsigned int m_num_n1, m_num_n2, m_num_n3;
    int m_adj_element[3];
    double m_phi0;
    double m_area;

    Node* m_node1;
    Node* m_node2;
    Node* m_node3;
    Edge*  m_edges[3];
    Hinge* m_hinges[3];
};

inline double Element::get_phi0() const { return m_phi0; }
inline double Element::get_area() const { return m_area; }

#endif //PLATES_SHELLS_ELEMENT_H
