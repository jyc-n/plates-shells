#ifndef PLATES_SHELLS_SIMULATION_H
#define PLATES_SHELLS_SIMULATION_H

#include <string>

class Parameters;
class Geometry;
class Boundary;
class PreProcessorImpl;
class SolverImpl;

class Simulation {

public:
    Simulation();
    ~Simulation();

    void pre_process(
                    const bool AR_FLAG, const bool K_FLAG,
                    int num_len, int num_wid,
                    double k_s, double k_sh, double k_b
                    );
    void solve();

private:
    Parameters* m_SimPar;
    Geometry*   m_SimGeo;
    Boundary*   m_SimBC;

    // Pointers to Implementations
    PreProcessorImpl* m_PreProcessor;
    SolverImpl*       m_SolverImpl;
};

#endif //PLATES_SHELLS_SIMULATION_H
