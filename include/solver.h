#ifndef PLATES_SHELLS_SOLVER_H
#define PLATES_SHELLS_SOLVER_H

#include "type_alias.h"

class Parameters;
class Node;
class Element;
class Geometry;
class Boundary;

/* PARDISO prototype. */
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                            double *, int    *,    int *, int *,   int *, int *,
                            int *, double *, double *, int *, double *);

const int SOLVER_TYPE = 1;      // 0 - Eigen CG solver, 1 - Pardiso (only works on linux)

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
    void findDofnew(const VectorN& rhs, SpMatrix& jacobian, VectorNodes& x_new);

    // solver functions
    void sparseSolver(const SpMatrix& A, const VectorN& rhs, VectorN& u);
    void pardisoInterface(SpMatrix& A, const VectorN& rhs, VectorN& u);

    // helper functions
    void findMappingVectors();
    void DEStretch(VectorN& dEdq, SparseEntries& entries_full);
    void DEShear  (VectorN& dEdq, SparseEntries& entries_full);
    void DEBend   (VectorN& dEdq, SparseEntries& entries_full);

    void analyticalStatic();
};

#endif //PLATES_SHELLS_SOLVER_H
