#include <iostream>
#include <cstdlib>
#include "simulation.h"

int main(int argc, char* argv[]) {
    try {
        bool AR_OVERRIDE = false, K_OVERRIDE = false;
        int num_len = 0, num_wid = 0;
        double k_s = 0.0, k_sh = 0.0, k_b = 0.0;
        switch (argc) {
            // every parameters read from input.txt
            case 1:
                break;

            // given #nodes along length and width
            case 3:
                num_len = atoi(argv[1]);
                num_wid = atoi(argv[2]);
                AR_OVERRIDE = true;
                break;

            // given #nodes along length and width, and k_s k_sh k_b
            case 6:
                num_len = atoi(argv[1]);
                num_wid = atoi(argv[2]);
                k_s     = atof(argv[3]);
                k_sh    = atof(argv[4]);
                k_b     = atof(argv[5]);
                AR_OVERRIDE = true;
                K_OVERRIDE = true;
                break;
        
            default:
                throw "wrong arguments!";
        }
        Simulation SimCase;
        SimCase.pre_process(AR_OVERRIDE, K_OVERRIDE, num_len, num_wid, k_s, k_sh, k_b);
        SimCase.solve();
    }
    catch (const char* msg) {
        std::cerr << msg << std::endl;
        exit(1);
    }
    return 0;
}