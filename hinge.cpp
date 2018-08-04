#include <iostream>
#include <cmath>
#include "node.h"
#include "element.h"
#include "hinge.h"

Hinge::Hinge()
 : m_el1(nullptr), m_el2(nullptr),
   m_node0(nullptr), m_node1(nullptr),
   m_node2(nullptr), m_node3(nullptr),
   m_n1n20(0), m_e0(0), VISITED(false)
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

void Hinge::find_n1n20() {
    m_n1n20 = ( m_el1->get_normal().norm() ) * ( m_el2->get_normal().norm() );
}

void Hinge::find_kappa0() {
    m_kappa0 = 2.0 * m_el1->get_normal().cross(m_el2->get_normal()) /
               ( m_n1n20 + m_el1->get_normal().dot(m_el2->get_normal()) );
}

void Hinge::find_e0() {
    Eigen::Vector3d m_eVec = m_node1->get_xyz() - m_node2->get_xyz();
    m_e0 = m_eVec.norm();
}

bool Hinge::check_visited() const {
    return VISITED;
}

double Hinge::get_area_sum() const {
    return m_el1->get_area() + m_el2->get_area();
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

double Hinge::get_n1n20() const {
    return m_n1n20;
}

Eigen::Vector3d Hinge::get_kappa0() const {
    return m_kappa0;
}

double Hinge::get_e0() const {
    return m_e0;
}