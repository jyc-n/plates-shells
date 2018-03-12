#ifndef PLATES_SHELLS_SOLVER_H
#define PLATES_SHELLS_SOLVER_H

#include "parameters.h"
#include "geometry.h"

void init_solver(Geometry& InitGeo, Parameters& SimPar);
void time_solver();

#endif //PLATES_SHELLS_SOLVER_H
