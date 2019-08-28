#include <iostream>
#include <cstdlib>
#include "arguments.h"
#include "simulation.h"

int main(int argc, char* argv[]) {
    try {
        Arguments t_args(argc, argv);

        Simulation SimCase(t_args.inputPath, t_args.outputPath);
        SimCase.pre_process(t_args);
        SimCase.solve();
    }
    catch (const char* msg) {
        std::cerr << msg << std::endl;
        exit(1);
    }
    return 0;
}