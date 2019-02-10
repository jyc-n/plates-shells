#include <iostream>
#include <cmath>
#include "element.h"
#include "node.h"
#include "edge.h"
#include "hinge.h"

Element::Element()
        : m_num_el(0),
          m_num_n1(0),
          m_num_n2(0),
          m_num_n3(0),
          m_node1(nullptr),
          m_node2(nullptr),
          m_node3(nullptr)
{
    for (int &i : m_adj_element)
        i = 0;
    for (int &i : m_edgeIndex)
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

    calculate_area();
    calculate_phi0();

    for (int &i : m_adj_element)
        i = 0;
    for (int &i : m_edgeIndex)
        i = 0;
    for (auto &i : m_edges)
        i = nullptr;
    for (auto &i : m_hinges)
        i = nullptr;
}

// NOTE: dynamically allocated Hinge and Edge classes are destroyed in Geometry destructor
Element::~Element()
{}

void Element::calculate_area() {
    Eigen::Vector3d dir12 = *((*m_node2).get_xyz()) - *((*m_node1).get_xyz());
    Eigen::Vector3d dir23 = *((*m_node3).get_xyz()) - *((*m_node2).get_xyz());
    m_area = 0.5 * (dir12.cross(dir23)).norm();
}

void Element::calculate_phi0() {
    Eigen::Vector3d dir12 = *((*m_node1).get_xyz()) - *((*m_node2).get_xyz());
    Eigen::Vector3d dir23 = *((*m_node3).get_xyz()) - *((*m_node2).get_xyz());
    dir12 = dir12 / dir12.norm();
    dir23 = dir23 / dir23.norm();
    m_phi0 = acos(dir12.dot(dir23));
}

void Element::find_edge_index() {
    int n1 = m_node1->get_num();
    int n2 = m_node2->get_num();
    int n3 = m_node3->get_num();
    for (int p = 0; p < 3; p++) {
        int i = 0, j = 0;
        switch (p) {
            case 0:
                i = std::max(n1, n2);
                j = std::min(n1, n2);
                break;
            case 1:
                i = std::max(n2, n3);
                j = std::min(n2, n3);
                break;
            case 2:
                i = std::max(n1, n3);
                j = std::min(n1, n3);
                break;
        }
        m_edgeIndex[p] = i*(i-1)/2 + j-1;
    }
}

void Element::find_adjacent_element(int index, int num_element) {
    // store the number of element that share the edge with given index
    m_adj_element[get_which_edge(index)] = num_element;
}

// construct the edge object and return the pointer to the edge
Edge* Element::build_edges(int index) {
    // find the local number of edge
    int local_number = get_which_edge(index);
    // build the edge with given index
    Node* n1 = nullptr;
    Node* n2 = nullptr;

    //  numbering convention of edges
    //            3
    //           /|
    //     e2  /  |
    //       /    | e1
    //     /______|
    //    1   e0  2

    switch (local_number) {
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

    m_edges[local_number] = new Edge(n1, n2);
//    std::cout << n1->get_num() << '\t' << n2->get_num() << '\n';
    return m_edges[local_number];
}

// construct the hinge object and return the pointer to the hinge
Hinge* Element::build_hinges(int index, Element* adj_element) {
    // find the local number of hinge (the overlapped edge)
    int local_number = get_which_edge(index);
    // build the hinge
    Node* n0 = nullptr;
    Node* n1 = nullptr;
    Node* n2 = nullptr;
    Node* n3 = nullptr;

    int n_remain_node = adj_element->get_remain_node_num(local_number);

    //  numbering convention of hinges
    //            2                  
    //           /|                        e4                       e3
    //     e1  /  |                    1_______ 3              2 ________ 1
    //       /    | e3                /|      /                 |      /|
    //   0 /______|             e3  /  |    /                   |  e0/  |  e4
    //    |   e0 /  1             /  e0|  / e2              e1  |  /    |
    //    |    /                /______|/                       |/______|
    // e2 |  /   e4           2   e1   0                       0    e2   3
    //    |/
    //    3
    //      case0                    case1                        case2
    // 
    // Edge orientation: e0,e1,e2 point away from node0
    //                      e3,e4 point away from node1

    switch (local_number) {
        case 0:
            n0 = m_node1;
            n1 = m_node2;
            n2 = m_node3;
            n3 = adj_element->get_node(n_remain_node);
            break;
        case 1:
            n0 = m_node2;
            n1 = m_node3;
            n2 = m_node1;
            n3 = adj_element->get_node(n_remain_node);
            break;
        case 2:
            n0 = m_node1;
            n1 = m_node3;
            n2 = adj_element->get_node(n_remain_node);
            n3 = m_node2;
            break;
    }

    m_hinges[local_number] = new Hinge(n0, n1, n2, n3, this, adj_element);
//    std::cout << m_hinges[local_number]->get_node(0)->get_num() << '\t'
//              << m_hinges[local_number]->get_node(1)->get_num() << '\t'
//              << m_hinges[local_number]->get_node(2)->get_num() << '\t'
//              << m_hinges[local_number]->get_node(3)->get_num() << '\n';
    return m_hinges[local_number];
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

int Element::get_edge_index(int num) const {
    switch (num) {
        case 0:
            return m_edgeIndex[0];
        case 1:
            return m_edgeIndex[1];
        case 2:
            return m_edgeIndex[2];
        default:
            std::cerr << "Local number of edge can only be 0, 1, 2" << std::endl;
            exit(1);
    }
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

int Element::get_which_edge(int index) const {
    // find the local number of edge with given index
    for (int i = 0; i < 3; i++) {
        if (index == m_edgeIndex[i])
            return i;
    }
    std::cerr << "Edge index " << index << " doesn't exists in element " << m_num_el << std::endl;
    exit(1);
}

// obtain the number of node that is NOT on the overlapped edge
int Element::get_remain_node_num(int num_edge) const {
    switch (num_edge) {
        case 0:
            return 3;
        case 1:
            return 1;
        case 2:
            return 2;
        default:
            std::cerr << "Local number of edge can only be 0, 1, 2" << std::endl;
            exit(1);
    }
}