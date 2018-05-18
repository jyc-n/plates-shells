#include "node.h"

Node::Node() {
    m_num = 0;
    m_x = 0.0;
    m_y = 0.0;
    m_z = 0.0;
}

Node::Node(const unsigned int& num, const double& xpos, const double& ypos, const double& zpos)
        : m_num(num), m_x(xpos), m_y(ypos), m_z(zpos)
{}

unsigned int Node::get_num() const {
    return m_num;
}

Eigen::Vector3d Node::get_xyz() const {
    Eigen::Vector3d result;
    result[0] = m_x;
    result[1] = m_y;
    result[2] = m_z;
    return result;
}