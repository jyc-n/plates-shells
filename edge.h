#ifndef PLATES_SHELLS_EDGE_H
#define PLATES_SHELLS_EDGE_H

class Node;

class Edge {
public:
    // TODO: write a different constructor, don't need the two "find" functions
    // constructor
    Edge(Node* n1, Node* n2);
    // destructor
    ~Edge();

    // accessor
    Node* get_node(int num) const;
    unsigned int get_node_num(const int num) const;
    double get_len0() const;

    double m_k;

private:
    Node* m_node1;
    Node* m_node2;

    double m_len0;
};

#endif //PLATES_SHELLS_EDGE_H
