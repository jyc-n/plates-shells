#ifndef PLATES_SHELLS_GEOMETRY_H
#define PLATES_SHELLS_GEOMETRY_H

#include "parameters.h"

class Node;
class Element;

class Geometry {

public:
    Geometry();

    // modifier
    void init_node_lst();
    void init_element_lst();

    // accessor
    Node& get_node_lst() const;
    Element& get_element_lst() const;

private:
    Node* m_node_lst;
    Element* m_element_lst;
};

#endif //PLATES_SHELLS_GEOMETRY_H
