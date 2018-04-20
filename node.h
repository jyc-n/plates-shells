#ifndef PLATES_SHELLS_NODE_H
#define PLATES_SHELLS_NODE_H

#include <Eigen/Dense>

class Node {
public:
    //constructor
    Node();
    Node(const double& xpos, const double& ypos, const double& zpos);

    // accessor
    Eigen::Vector3d get_xyz() const;

private:
    double m_x, m_y, m_z;
};

#endif //PLATES_SHELLS_NODE_H
