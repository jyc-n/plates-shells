#include <iostream>
#include <cmath>
#include "node.h"
#include "element.h"
#include "hinge.h"

Hinge::Hinge()
 : m_el1(nullptr), m_el2(nullptr),
   m_node0(nullptr), m_node1(nullptr),
   m_node2(nullptr), m_node3(nullptr),
   VISITED(false)
{}

void Hinge::mark_visited() {
    VISITED = true;
}

void Hinge::reset_visited() {
    VISITED = false;
}

void Hinge::find_element(Element& el1, Element& el2) {
    m_el1 = &el1;
    m_el2 = &el2;
}

void Hinge::find_node(Node& n0, Node& n1, Node& n2, Node& n3) {
    m_node0 = &n0;
    m_node1 = &n1;
    m_node2 = &n2;
    m_node3 = &n3;
}

void Hinge::find_originVal(){
    Eigen::Vector3d ve0 = m_node1->get_xyz() - m_node0->get_xyz();
    double e0 = ve0.norm();
    double a0 = m_el1->get_area() + m_el2->get_area();
    m_const = 3.0 * pow(e0, 2) / a0;

    //TODO: calculate psi0
    m_psi0 = 0;
}

bool Hinge::check_visited() const {
    return VISITED;
}

Node* Hinge::get_node(int num) const {
    switch (num) {
        case 0:
            return m_node0;
        case 1:
            return m_node1;
        case 2:
            return m_node2;
        case 3:
            return m_node3;
        default:
            std::cerr << "Local number of node can only be 0, 1, 2, 3" << std::endl;
            exit(1);
    }
}

unsigned int Hinge::get_node_num(const int num) const {
    switch (num) {
        case 0:
            return m_node0->get_num();
        case 1:
            return m_node1->get_num();
        case 2:
            return m_node2->get_num();
        case 3:
            return m_node3->get_num();
        default:
            std::cerr << "Local number of node can only be 0, 1, 2, 3" << std::endl;
            exit(1);
    }
}

