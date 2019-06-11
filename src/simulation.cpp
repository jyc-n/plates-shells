#include <iostream>
#include "simulation.h"
#include "arguments.h"
#include "parameters.h"
#include "geometry.h"
#include "loadbc.h"
#include "pre_processor.h"
#include "solver.h"
#include "utilities.h"

Simulation::Simulation(const std::string& t_input, const std::string& t_output) {
    m_SimPar = new Parameters(t_input, t_output);
    m_SimGeo = new Geometry(m_SimPar);
    m_SimBC = new Boundary(m_SimPar, m_SimGeo);

    m_PreProcessor = new PreProcessorImpl(m_SimPar, m_SimGeo);
    m_SolverImpl = new SolverImpl(m_SimPar, m_SimGeo, m_SimBC);
}

Simulation::~Simulation() {
    delete m_SimPar;
    delete m_SimGeo;
    delete m_SimBC;
    delete m_PreProcessor;
    delete m_SolverImpl;
}

void Simulation::pre_process(const Arguments& t_args) {
    m_PreProcessor->PreProcess(t_args);
    m_SimBC->initBC();
    m_SolverImpl->initSolver();
}

void Simulation::solve() {
    if (SOLVER_TYPE == 0)
        std::cout << "Eigen CG solver will be used" << std::endl;
    else if (SOLVER_TYPE == 1)
        std::cout << "Pardiso solver will be used" << std::endl;

    // timer for entire simulation
    Timer t_all(true);
    
    if (m_SimPar->solver_op()) {
        std::cout << "Dynamic simulation starts\n" << std::endl;
        m_SolverImpl->dynamic();
    }
    else {
        std::cout << "Static simulation starts\n" << std::endl;
        m_SolverImpl->statics();
    }

    std::cout << "---------------------------" << std::endl;
    std::cout << "Simulation completed" << std::endl;
    std::cout << "Total time used " << t_all.elapsed(true) << " seconds" <<  std::endl;
}
