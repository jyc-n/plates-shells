#include "pre_processor.h"
#include "solver.h"

int main() {
    try {
        pre_processor();
        time_solver();
    }
    catch (const char* msg) {
        std::cerr << msg << std::endl;
        exit(1);
    }
    return 0;
}