#ifndef PLATES_SHELLS_HINGE_H
#define PLATES_SHELLS_HINGE_H

#include <Eigen/Dense>

class Node;
class Element;

class Hinge {
public:

    // constructor
    Hinge();

    // modifiers
    void mark_visited();
    void reset_visited();
    void find_element(Element& el1, Element& el2);
    void find_node(Node& n0, Node& n1, Node& n2, Node& n3);
    void find_n1n20();
    void find_kappa0();
    void find_e0();

    // accessor
    bool check_visited() const;
    double get_area_sum() const;
    Node* get_node(int num) const;
    unsigned int get_node_num(const int num) const;
    double get_n1n20() const;
    Eigen::Vector3d get_kappa0() const;
    double get_e0() const;

private:
    Element* m_el1;
    Element* m_el2;

    Node* m_node0;
    Node* m_node1;
    Node* m_node2;
    Node* m_node3;

    Eigen::Vector3d m_kappa0;
    double m_n1n20;
    double m_e0;
    bool VISITED;
};

#endif //PLATES_SHELLS_HINGE_H
