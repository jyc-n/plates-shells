#include <iostream>
#include <cmath>
#include "geometry.h"
#include "element.h"
#include "node.h"
#include "edge.h"
#include "hinge.h"

Element::Element()
        : m_num_el(0),
          m_num_n1(0),
          m_num_n2(0),
          m_num_n3(0)
{
    for (int &i : m_adj_element)
        i = 0;
    for (auto &i : m_edges)
        i = nullptr;
    for (auto &i : m_hinges)
        i = nullptr;
}

Element::Element(const unsigned int num_el, Node* n1, Node* n2, Node* n3)
        : m_num_el(num_el),
          m_node1(n1),
          m_node2(n2),
          m_node3(n3)
{
    m_num_n1 = (*m_node1).get_num();
    m_num_n2 = (*m_node2).get_num();
    m_num_n3 = (*m_node3).get_num();

    calculate_vec_edge();
    calculate_len_edge();
    calculate_area();
    calculate_normal();

    for (int &i : m_adj_element)
        i = 0;
    for (auto &i : m_edges)
        i = nullptr;
    for (auto &i : m_hinges)
        i = nullptr;
}

// NOTE: dynamically allocated Hinge and Edge classes are destroyed in Geometry destructor
Element::~Element()
{}

void Element::calculate_vec_edge() {
    m_dir12 = (*m_node2).get_xyz() - (*m_node1).get_xyz();
    m_dir23 = (*m_node3).get_xyz() - (*m_node2).get_xyz();
    m_dir13 = (*m_node3).get_xyz() - (*m_node1).get_xyz();
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

void Element::find_nearby_element(Geometry& Geo) {
    for (int i = 0; i < Geo.nel(); i++) {    // looking at itself in the element list
        if (m_num_el == i+1) {
            continue;
        }

        bool N1_OVERLAP = false;
        bool N2_OVERLAP = false;
        bool N3_OVERLAP = false;

        for (int j = 0; j < Geo.ndof(); j++) {
            if (m_num_n1 == Geo.m_conn(i,0) || m_num_n1 == Geo.m_conn(i,1) || m_num_n1 == Geo.m_conn(i,2))
                N1_OVERLAP = true;
            if (m_num_n2 == Geo.m_conn(i,0) || m_num_n2 == Geo.m_conn(i,1) || m_num_n2 == Geo.m_conn(i,2))
                N2_OVERLAP = true;
            if (m_num_n3 == Geo.m_conn(i,0) || m_num_n3 == Geo.m_conn(i,1) || m_num_n3 == Geo.m_conn(i,2))
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

void Element::find_edges(std::vector<Element>& l_element) {
    for (int i = 0; i < 3; i++) {

        // if the edge is a hinge (shared by two elements)
        if (this->is_hinge(i+1)) {

            Element* adj_element = &l_element[m_adj_element[i]-1];
            int n_edge = adj_element->get_overlapped_edge_num(m_num_el);

            // check if the hinge has already been set up
            if (adj_element->get_edge(n_edge) != nullptr) {
                m_edges[i] = adj_element->get_edge(n_edge);
                continue;
            }
        }
        // if this edge is not a hinge OR it's a hinge but the edge has not been set up yet
        m_edges[i] = new Edge;
        Node* n1 = nullptr;
        Node* n2 = nullptr;

        switch (i) {
            case 0:
                n1 = m_node1;
                n2 = m_node2;
                break;
            case 1:
                n1 = m_node2;
                n2 = m_node3;
                break;
            case 2:
                n1 = m_node3;
                n2 = m_node1;
                break;
        }
        m_edges[i]->find_node(*n1, *n2);
        m_edges[i]->find_len0();
    }
}

void Element::find_hinges(std::vector<Element>& l_element) {
    for (int i = 0; i < 3; i++) {

        if (m_adj_element[i] != 0) {

            Element* adj_element = &l_element[m_adj_element[i]-1];

            int n_edge = adj_element->get_overlapped_edge_num(m_num_el);

            // if the hinge has already set up
            if (adj_element->get_hinge(n_edge) != nullptr)
                m_hinges[i] = adj_element->get_hinge(n_edge);

            else {
                m_hinges[i] = new Hinge;
                Node* n0 = nullptr;
                Node* n1 = nullptr;
                Node* n2 = nullptr;

                m_hinges[i]->find_element(*this, *adj_element);

                int n_remain_node = adj_element->get_remain_node_num(n_edge);
                Node* n3 = adj_element->get_node(n_remain_node);

                switch (i) {
                    case 0:
                        n0 = m_node2;
                        n1 = m_node1;
                        n2 = m_node3;
                        break;
                    case 1:
                        n0 = m_node2;
                        n1 = m_node3;
                        n2 = m_node1;
                        break;
                    case 2:
                        n0 = m_node1;
                        n1 = m_node3;
                        n2 = m_node2;
                        break;
                }
                m_hinges[i]->find_node(*n0, *n1, *n2, *n3);
                m_hinges[i]->find_n1n20();
                m_hinges[i]->find_kappa0();
                m_hinges[i]->find_e0();
            }
        }
    }
}

unsigned int Element::get_node_num(const int num) const {
    switch (num) {
        case 1:
            return m_num_n1;
        case 2:
            return m_num_n2;
        case 3:
            return m_num_n3;
        default:
            std::cerr << "Local number of node can only be 1, 2, 3" << std::endl;
            exit(1);
    }
}

double Element::get_area() const {
    return m_area;
}

double Element::get_len_edge(const int num) const {
    if (!(num == 1 || num == 2 || num == 3)){
        std::cerr << "Local number of edge can only be 1, 2, 3" << std::endl;
        exit(1);
    }
    return m_len_edge[num-1];
}

Node* Element::get_node(const int num) const {
    switch (num) {
        case 1:
            return m_node1;
        case 2:
            return m_node2;
        case 3:
            return m_node3;
        default:
            std::cerr << "Local number of node can only be 1, 2, 3" << std::endl;
            exit(1);
    }
}

int Element::get_nearby_element(int num) const {
    // indices of the edge should be 1, 2, 3
    num -= 1;
    if (num < 0 || num > 2) {
        std::cerr << "Local number of edge can only be 1, 2, 3" << std::endl;
        exit(1);
    }
    return m_adj_element[num];
}

// obtain the number of overlapped edge of the adjacent element
int Element::get_overlapped_edge_num(int num_el) const {
    for (int edge = 0; edge < 3; edge++) {
        if (m_adj_element[edge] == num_el)
            return edge+1;
    }
    return 0;
}

// obtain the number of node that is NOT on the overlapped edge
int Element::get_remain_node_num(int num_edge) const {
    switch (num_edge) {
        case 1:
            return 3;
        case 2:
            return 1;
        case 3:
            return 2;
        default:
            std::cerr << "Local number of edge can only be 1, 2, 3" << std::endl;
            exit(1);
    }
}

// check if an edge is a hinge
bool Element::is_hinge(int num_edge) const {
    if (num_edge < 1 || num_edge > 3) {
        std::cerr << "Local number of edge can only be 1, 2, 3" << std::endl;
        exit(1);
    }
    return (m_adj_element[num_edge - 1] != 0);
}

Edge* Element::get_edge(int num_edge) const {
    if (num_edge < 1 || num_edge > 3) {
        std::cerr << "Local number of edge can only be 1, 2, 3" << std::endl;
        exit(1);
    }
    return m_edges[num_edge-1];
}

// get the pointer to the hinge
Hinge* Element::get_hinge(int num_edge) const {

    // if the edge is a hinge
    if (this->is_hinge(num_edge))
        return m_hinges[num_edge-1];
    return nullptr;
}

Eigen::Vector3d Element::get_vec_edge(const int num) const {
    switch (num) {
        case 1:
            return m_dir12;
        case 2:
            return m_dir23;
        case 3:
            return m_dir13;
        default:
            std::cerr << "Local number of edge can only be 1, 2, 3" << std::endl;
            exit(1);
    }
}

Eigen::Vector3d Element::get_normal() const {
    return m_normal;
}