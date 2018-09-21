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
    void find_originVal();

    // accessor
    bool check_visited() const;
    Node* get_node(int num) const;
    unsigned int get_node_num(const int num) const;

    double m_psi0;
    double m_const;
    double m_k;

private:
    Element* m_el1;
    Element* m_el2;

    Node* m_node0;
    Node* m_node1;
    Node* m_node2;
    Node* m_node3;

    bool VISITED;
};

#endif //PLATES_SHELLS_HINGE_H
