#ifndef PLATES_SHELLS_EDGE_H
#define PLATES_SHELLS_EDGE_H

#include <Eigen/Dense>

class Node;

class Edge {
public:

    // constructor
    Edge();

    // modifiers
    void mark_visited();
    void reset_visited();
    void find_node(Node& n1, Node& n2);
    void find_originVal();

    // accessor
    bool check_visited() const;
    Node* get_node(int num) const;
    unsigned int get_node_num(const int num) const;
    double get_len0() const;

    double m_k;

private:
    Node* m_node1;
    Node* m_node2;
    Eigen::Vector3d m_eVec;

    double m_len0;
    bool VISITED;
};

#endif //PLATES_SHELLS_EDGE_H
