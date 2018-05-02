#include <iostream>
#include <vector>
#include "pre_processor.h"
#include "solver.h"

int main() {
    try {
        Parameters SimParams;
        pre_processor(SimParams);
        //fake_solver(SimParams);
        time_solver(SimParams);
    }
    catch (const char* msg) {
        std::cerr << msg << std::endl;
        exit(1);
    }
    return 0;
}