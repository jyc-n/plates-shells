#ifndef PLATES_SHELLS_ELEMENT_H
#define PLATES_SHELLS_ELEMENT_H

#include <vector>
#include <Eigen/Dense>
#include "node.h"

class Element {
public:

    // constructor
    Element(const unsigned int& num_el,
            const unsigned int& num_n1,
            const unsigned int& num_n2,
            const unsigned int& num_n3);

    // modifier
    void set_node(std::vector<Node>& vec);
    void calculate_dir();
    void calculate_angle();
    void calculate_normal();

    // accessor
    Node& get_node(const unsigned int& num);
    double get_angle(const unsigned int& num);
    double get_normal(const unsigned int& num);

private:
    unsigned int m_num_el;
    unsigned int m_num_n1, m_num_n2, m_num_n3;
    double m_angle1, m_angle2, m_angle3;
    Eigen::Vector3d m_normal, m_dir12, m_dir23, m_dir13;

    Node node1;
    Node node2;
    Node node3;
};

#endif //PLATES_SHELLS_ELEMENT_H
