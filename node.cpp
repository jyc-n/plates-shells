#include "node.h"

Node::Node() {
    m_x = 0.0;
    m_y = 0.0;
    m_z = 0.0;
}

Node::Node(const double& xpos,
           const double& ypos,
           const double& zpos)
        : m_x(xpos), m_y(ypos), m_z(zpos)
{}

Eigen::Vector3d Node::get_xyz() const {
    Eigen::Vector3d result;
    result[0] = m_x;
    result[1] = m_y;
    result[2] = m_z;
    return result;
}

