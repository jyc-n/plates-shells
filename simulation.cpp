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

void Simulation::pre_process() {
    m_PreProcessor->PreProcess();
    m_SimBC->initBC();
}

void Simulation::solve() {
    m_SolverImpl->Solve();
}
