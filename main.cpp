#include <iostream>
#include <vector>
#include "pre_processor.h"
#include "solver.h"
#include "element.h"

int main() {
    try {
        Parameters SimParams;
        Geometry SimGeo;
        std::vector<Element> initialGeo = pre_processor(SimGeo, SimParams);
        fake_solver(SimGeo, SimParams);
    }
    catch (const char* msg) {
        std::cerr << msg << std::endl;
        exit(1);
    }
    return 0;
}