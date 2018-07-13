#ifndef PLATES_SHELLS_NODE_H
#define PLATES_SHELLS_NODE_H

#include <Eigen/Dense>

class Node {
public:
    //constructor
    Node();
    Node(const unsigned int num, const double xpos, const double ypos, const double zpos);

    // modifier
    void set_xyz(double x, double y, double z);

    // accessor
    unsigned int get_num() const;
    Eigen::Vector3d get_xyz() const;

private:
    unsigned int m_num;
    double m_x, m_y, m_z;
};

#endif //PLATES_SHELLS_NODE_H
