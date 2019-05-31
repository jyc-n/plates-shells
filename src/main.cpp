#include <iostream>
#include <cstdlib>
#include "arguments.h"
#include "simulation.h"

int main(int argc, char* argv[]) {
    try {
        Arguments t_args(argc);
        switch (argc) {
            // every parameters read from input.txt
            case DEFAULT:
                break;

            // given #nodes along length and width
            case NUMS:
                t_args.getArgs(atoi(argv[1]), atoi(argv[2]));
                break;

            // given #nodes along length and width, and length,width
            case NUMS_DIMS:
                t_args.getArgs(atoi(argv[1]), atoi(argv[2]), atof(argv[3]), atof(argv[4]));
                break;
        
            default:
                throw "wrong arguments!";
        }

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