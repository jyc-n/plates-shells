#ifndef PLATES_SHELLS_SOLVER_H
#define PLATES_SHELLS_SOLVER_H

#include <fstream>
#include "parameters.h"

class Node;
class Element;

// different solver
void fake_solver(Parameters& SimPar);
void time_solver(Parameters& SimPar);

// other functions
void prepare_solver(Parameters& SimPar, Node* l_node, Element* l_element);
//void write_output(std::string& filename, Eigen::MatrixXd& mat);

#endif //PLATES_SHELLS_SOLVER_H
