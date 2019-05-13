#ifndef PLATES_SHELLS_SIMULATION_H
#define PLATES_SHELLS_SIMULATION_H

#include <string>

class Arguments;
class Parameters;
class Geometry;
class Boundary;
class PreProcessorImpl;
class SolverImpl;

class Simulation {

public:
    Simulation(const std::string& t_input, const std::string& t_output);
    ~Simulation();

    void pre_process(Arguments t_args);
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
