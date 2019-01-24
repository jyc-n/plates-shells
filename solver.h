#ifndef PLATES_SHELLS_SOLVER_H
#define PLATES_SHELLS_SOLVER_H

#include "type_alias.h"

class Parameters;
class Node;
class Element;
class Geometry;
class Boundary;

class SolverImpl {

public:
    SolverImpl(Parameters* SimPar, Geometry* SimGeo, Boundary* SimBC);

    void initSolver();

    // main solver function
    bool step(const int ist, VectorNodes& x, VectorNodes& x_new, VectorNodes& vel);
    void dynamic();
    void quasistatic();
    void writeToFiles(const int ist);

private:

    // pointers
    Parameters* m_SimPar;
    Geometry*   m_SimGeo;
    Boundary*   m_SimBC;

    // member variables
    unsigned int m_numTotal;
    unsigned int m_numDirichlet;
    unsigned int m_numNeumann;
    double m_tol;
    double m_incRatio;

    std::vector<int> m_fullToDofs;
    std::vector<int> m_dofsToFull;
    
    // solver flags
    bool DYNAMIC_SOLVER;        // true - dynamic solver, false - static solver
    bool WRITE_OUTPUT;          // true - write output, false - no output, configured in input.txt
    bool INFO_STYLE;            // true - print info on console, false - progress bar + log file

    // subroutine
    void findDEnergy(VectorN& dEdq, SparseEntries& entries_full);
    void findResidual(const VectorNodes& vel, const VectorNodes& x, const VectorNodes& x_new, const VectorN& dEdq, VectorN& rhs);   // dynamic
    void findJacobian(SparseEntries& entries_full, SpMatrix& jacobian);
    void findDofnew(const VectorN& rhs, const SpMatrix& jacobian, VectorNodes& x_new);

    // solver functions
    void sparseSolver(const SpMatrix& A, const VectorN& b, VectorN& x);

    // helper functions
    void findMappingVectors();
    void DEStretch(VectorN& dEdq, SparseEntries& entries_full);
    void DEShear  (VectorN& dEdq, SparseEntries& entries_full);
    void DEBend   (VectorN& dEdq, SparseEntries& entries_full);
    void calcViscous(const Eigen::VectorXd& qn, const Eigen::VectorXd& qnew, Eigen::VectorXd& dEdq, Eigen::MatrixXd& ddEddq);

    void analyticalStatic();
};

#endif //PLATES_SHELLS_SOLVER_H
