#ifndef PLATES_SHELLS_SOLVER_H
#define PLATES_SHELLS_SOLVER_H

#include <vector>
#include <Eigen/Dense>

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
    void Solve();

private:

    // subroutine
    void calcDEnergy(Eigen::VectorXd& dEdq, Eigen::MatrixXd& ddEddq);
    void updateResidual(Eigen::VectorXd& qn, Eigen::VectorXd& qnew, Eigen::VectorXd& vel, Eigen::VectorXd& dEdq, Eigen::VectorXd& res_f);
    void updateJacobian(Eigen::MatrixXd& ddEddq, Eigen::MatrixXd& mat_j);
    Eigen::VectorXd calcVel(double dt, const Eigen::VectorXd& qcurr, const Eigen::VectorXd& qnew);
    Eigen::VectorXd calcDofnew(Eigen::VectorXd& qn, const Eigen::VectorXd& temp_f, const Eigen::MatrixXd& temp_j);

    // solver functions
    void denseSolver(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x);
    void sparseSolver(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x);

    // helper functions
    Eigen::VectorXd unconsVec(const Eigen::VectorXd& vec);
    Eigen::MatrixXd unconsMat(const Eigen::MatrixXd& mat);
    void updateNodes(Eigen::VectorXd& qnew);

    void calcStretch(Eigen::VectorXd& dEdq, Eigen::MatrixXd& ddEddq);
    void calcShear  (Eigen::VectorXd& dEdq, Eigen::MatrixXd& ddEddq);
    void calcBend   (Eigen::VectorXd& dEdq, Eigen::MatrixXd& ddEddq);
    void calcViscous(const Eigen::VectorXd& qn, const Eigen::VectorXd& qnew, Eigen::VectorXd& dEdq, Eigen::MatrixXd& ddEddq);

    // pointers
    Parameters* m_SimPar;
    Geometry*   m_SimGeo;
    Boundary*   m_SimBC;

    // member variables
    Eigen::VectorXd m_dEdq;
    Eigen::MatrixXd m_ddEddq;
    Eigen::VectorXd m_residual;
    Eigen::MatrixXd m_jacobian;
};

#endif //PLATES_SHELLS_SOLVER_H
