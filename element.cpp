#include <iostream>
#include <cmath>
#include "geometry.h"
#include "element.h"

Element::Element()
        : m_num_el(0),
          m_num_n1(0),
          m_num_n2(0),
          m_num_n3(0)
{
    for (int &i : m_adj_element)
        i = 0;
}

Element::Element(const unsigned int& num_el, Node* n1, Node* n2, Node* n3)
        : m_num_el(num_el),
          node1(n1),
          node2(n2),
          node3(n3)
{
    m_num_n1 = (*node1).get_num();
    m_num_n2 = (*node2).get_num();
    m_num_n3 = (*node3).get_num();

    calculate_vec_edge();
    calculate_len_edge();
    calculate_area();
    calculate_normal();

    for (int &i : m_adj_element)
        i = 0;
}

void Element::calculate_vec_edge() {
    m_dir12 = (*node2).get_xyz() - (*node1).get_xyz();
    m_dir23 = (*node3).get_xyz() - (*node2).get_xyz();
    m_dir13 = (*node3).get_xyz() - (*node1).get_xyz();
}

void Element::calculate_len_edge() {
    m_len_edge[0] = m_dir12.norm();
    m_len_edge[1] = m_dir23.norm();
    m_len_edge[2] = m_dir13.norm();
}

void Element::calculate_area() {
    m_area = 0.5 * (m_dir12.cross(m_dir23)).norm();
}

void Element::calculate_normal() {
    m_normal = m_dir12.cross(m_dir23);
}

void Element::find_nearby_element(const Parameters& Params) {
    for (int i = 0; i < Params.nel(); i++) {    // looking at itself in the element list
        if (m_num_el == i+1) {
            continue;
        }

        bool N1_OVERLAP = false;
        bool N2_OVERLAP = false;
        bool N3_OVERLAP = false;

        for (int j = 0; j < Params.ndof(); j++) {
            if (m_num_n1 == m_conn(i,0) || m_num_n1 == m_conn(i,1) || m_num_n1 == m_conn(i,2))
                N1_OVERLAP = true;
            if (m_num_n2 == m_conn(i,0) || m_num_n2 == m_conn(i,1) || m_num_n2 == m_conn(i,2))
                N2_OVERLAP = true;
            if (m_num_n3 == m_conn(i,0) || m_num_n3 == m_conn(i,1) || m_num_n3 == m_conn(i,2))
                N3_OVERLAP = true;
        }

        if (N1_OVERLAP && N2_OVERLAP)
            m_adj_element[0] = i+1;         // element number is the position in the array +1
        if (N2_OVERLAP && N3_OVERLAP)
            m_adj_element[1] = i+1;
        if (N1_OVERLAP && N3_OVERLAP)
            m_adj_element[2] = i+1;
    }
}

unsigned int Element::get_node_num(const int& num) const {
    if (!(num == 1 || num == 2 || num == 3))
        std::cerr << "Local number of node can only be 1, 2, 3" << std::endl;
    switch (num) {
        case 1:
            return m_num_n1;
        case 2:
            return m_num_n2;
        case 3:
            return m_num_n3;
    }
}

double Element::get_area() const {
    return m_area;
}

double Element::get_len_edge(const int& num) const {
    if (!(num == 1 || num == 2 || num == 3))
        std::cerr << "Local number of edge can only be 1, 2, 3" << std::endl;
    return m_len_edge[num-1];
}

Node* Element::get_node(const int& num) const {
    if (!(num == 1 || num == 2 || num == 3))
        std::cerr << "Local number of node can only be 1, 2, 3" << std::endl;
    switch (num) {
        case 1:
            return node1;
        case 2:
            return node2;
        case 3:
            return node3;
    }
}

Eigen::Vector3d Element::get_vec_edge(const int& num) const {
    if (!(num == 1 || num == 2 || num == 3))
        std::cerr << "Local number of edge can only be 1, 2, 3" << std::endl;
    switch (num) {
        case 1:
            return m_dir12;
        case 2:
            return m_dir23;
        case 3:
            return m_dir13;
    }
}

Eigen::Vector3d Element::get_normal() const {
    return m_normal;
}