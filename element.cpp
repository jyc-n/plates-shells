#include <cmath>
#include "geometry.h"
#include "element.h"

Element::Element()
        : m_num_el(0),
          m_num_n1(0),
          m_num_n2(0),
          m_num_n3(0)
{}

Element::Element(const unsigned int& num_el,
                 const unsigned int& num_n1,
                 const unsigned int& num_n2,
                 const unsigned int& num_n3)
        : m_num_el(num_el),
          m_num_n1(num_n1),
          m_num_n2(num_n2),
          m_num_n3(num_n3)
{}

void Element::set_node(Node* n1, Node* n2, Node* n3) {
    node1 = n1;
    node2 = n2;
    node3 = n3;
}

void Element::calculate_dir() {
    m_dir12 = (*node2).get_xyz() - (*node1).get_xyz();
    m_dir23 = (*node3).get_xyz() - (*node2).get_xyz();
    m_dir13 = (*node3).get_xyz() - (*node1).get_xyz();

    m_dir12 = m_dir12 / m_dir12.norm();
    m_dir23 = m_dir23 / m_dir23.norm();
    m_dir13 = m_dir13 / m_dir13.norm();
}

void Element::calculate_angle() {
    m_angle1 = acos(m_dir12.dot(m_dir13) / (m_dir12.norm() * m_dir13.norm()));
    m_angle2 = acos(m_dir12.dot(m_dir23) / (m_dir12.norm() * m_dir23.norm()));
    m_angle3 = acos(m_dir13.dot(m_dir23) / (m_dir13.norm() * m_dir23.norm()));
}

void Element::calculate_normal() {
    m_normal = m_dir12.cross(m_dir23);
}

void Element::find_nearby_element(const Parameters& Params) {
    for (int i = 0; i < Params.nel(); i++) {
        if (m_num_el == i+1) {
            continue;
        }

        if (m_num_n1 == m_conn(i,0) || m_conn(i,1) || m_conn(i,2) ) {
            if (m_num_n2 == m_conn(i,0) || m_conn(i,1) || m_conn(i,2)) {
                m_adj_element[0] = i;
            }
            else if (m_num_n3 == m_conn(i,0) || m_conn(i,1) || m_conn(i,2)) {
                m_adj_element[2] = i;
            }

        }
        else if (m_num_n2 == m_conn(i,0) || m_conn(i,1) || m_conn(i,2) ) {
            if (m_num_n1 == m_conn(i,0) || m_conn(i,1) || m_conn(i,2)) {
                m_adj_element[0] = i;
            }
            else if (m_num_n3 == m_conn(i,0) || m_conn(i,1) || m_conn(i,2)) {
                m_adj_element[1] = i;
            }
        }
        else if (m_num_n3 == m_conn(i,0) || m_conn(i,1) || m_conn(i,2) ) {
            if (m_num_n1 == m_conn(i,0) || m_conn(i,1) || m_conn(i,2)) {
                m_adj_element[2] = i;
            }
            else if (m_num_n2 == m_conn(i,0) || m_conn(i,1) || m_conn(i,2)) {
                m_adj_element[1] = i;
            }
        }
    }
}
