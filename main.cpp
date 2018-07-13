#include <iostream>
#include "simulation.h"

int main() {
    try {
        Simulation SimCase;
        SimCase.pre_process();
        SimCase.solve();
    }
    catch (const char* msg) {
        std::cerr << msg << std::endl;
        exit(1);
    }
    return 0;
}