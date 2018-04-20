#ifndef PLATES_SHELLS_SOLVER_H
#define PLATES_SHELLS_SOLVER_H

#include <fstream>
#include "parameters.h"

class Element;

void fake_solver(Parameters& SimPar);
void time_solver();

//void write_output(std::string& filename, Eigen::MatrixXd& mat);
#endif //PLATES_SHELLS_SOLVER_H
