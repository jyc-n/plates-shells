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

    void pre_process();
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
