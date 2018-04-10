#include "node.h"

Node::Node() {
    x_ = 0.0;
    y_ = 0.0;
    z_ = 0.0;
}

Node::Node(const double& xpos,
           const double& ypos,
           const double& zpos)
        : x_(xpos), y_(ypos), z_(zpos)
{}

Eigen::Vector3d Node::get_xyz() const {
    Eigen::Vector3d result;
    result[0] = x_;
    result[1] = y_;
    result[2] = z_;
    return result;
}

