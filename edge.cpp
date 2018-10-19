#include <iostream>
#include "node.h"
#include "edge.h"

Edge::Edge()
  : m_node1(nullptr), m_node2(nullptr),
    m_len0(0), VISITED(false)
{}

void Edge::mark_visited() {
    VISITED = true;
}

void Edge::reset_visited() {
    VISITED = false;
}

void Edge::find_node(Node &n1, Node &n2) {
    m_node1 = &n1;
    m_node2 = &n2;
}

void Edge::find_originVal() {
    m_eVec = m_node1->get_xyz() - m_node2->get_xyz();
    m_len0 = m_eVec.norm();
}

bool Edge::check_visited() const {
    return VISITED;
}

Node* Edge::get_node(int num) const {
    switch (num) {
        case 1:
            return m_node1;
        case 2:
            return m_node2;
        default:
            std::cerr << "Local number of node can only be 1, 2" << std::endl;
            exit(1);
    }
}

unsigned int Edge::get_node_num(const int num) const {
    switch (num) {
        case 1:
            return m_node1->get_num();
        case 2:
            return m_node2->get_num();
        default:
            std::cerr << "Local number of node can only be 1, 2" << std::endl;
            exit(1);
    }
}

double Edge::get_len0() const {
    return m_len0;
}