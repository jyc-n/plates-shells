#include "simulation.h"
#include "parameters.h"
#include "geometry.h"
#include "loadbc.h"
#include "pre_processor.h"
#include "solver.h"

Simulation::Simulation() {
    m_SimPar = new Parameters;
    m_SimGeo = new Geometry(m_SimPar);
    m_SimPar->link_geo(m_SimGeo);
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

void Simulation::pre_process(const bool AR_FLAG, const bool K_FLAG,
                             int num_len, int num_wid,
                             double k_s, double k_sh, double k_b) {
    m_PreProcessor->PreProcess(AR_FLAG, K_FLAG, num_len, num_wid, k_s, k_sh, k_b);
    m_SimBC->initBC();
    m_SolverImpl->initSolver();
}

void Simulation::solve() {
    if (m_SimPar->solver_op())
        m_SolverImpl->dynamic();
    else
        m_SolverImpl->quasistatic();
}
