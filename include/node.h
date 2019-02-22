#ifndef PLATES_SHELLS_NODE_H
#define PLATES_SHELLS_NODE_H

#include <Eigen/Dense>

class Node {
public:
    //constructor
    Node()
        : m_num(0), m_xyz(nullptr)
    {}
    Node(const unsigned int num, Eigen::Vector3d* xyz)
        : m_num(num), m_xyz(xyz)
    {}

    // accessor
    unsigned int get_num() const { return m_num; }
    Eigen::Vector3d* get_xyz() { return m_xyz; }

private:
    unsigned int m_num;
    Eigen::Vector3d* m_xyz;
};

#endif //PLATES_SHELLS_NODE_H
