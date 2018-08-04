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

    // subroutine
    void calcDEnergy(Eigen::VectorXd& dEdq, Eigen::MatrixXd& ddEddq);
    void updateResidual(Eigen::VectorXd& qn, Eigen::VectorXd& qnew, Eigen::VectorXd& vel, Eigen::VectorXd& dEdq, Eigen::VectorXd& res_f);
    void updateJacobian(Eigen::MatrixXd& ddEddq, Eigen::MatrixXd& mat_j);
    Eigen::VectorXd calcVel(double dt, const Eigen::VectorXd& qcurr, const Eigen::VectorXd& qnew);
    Eigen::VectorXd calcDofnew(Eigen::VectorXd& qn, const Eigen::VectorXd& temp_f, const Eigen::MatrixXd& temp_j);

    // helper functions
    Eigen::VectorXd unconsVec(const Eigen::VectorXd& vec);
    Eigen::MatrixXd unconsMat(const Eigen::MatrixXd& mat);
    void updateNodes(Eigen::VectorXd& qnew);

    void calcStretch(Eigen::VectorXd& dEdq, Eigen::MatrixXd& ddEddq);
    void calcShear  (Eigen::VectorXd& dEdq, Eigen::MatrixXd& ddEddq);
    void calcBend   (Eigen::VectorXd& dEdq, Eigen::MatrixXd& ddEddq);

private:
    Parameters* m_SimPar;
    Geometry*   m_SimGeo;
    Boundary*   m_SimBC;

    // member variables
    Eigen::VectorXd m_dEdq;
    Eigen::MatrixXd m_ddEddq;
    Eigen::VectorXd m_residual;
    Eigen::MatrixXd m_jacobian;
};

// other non-member helper functions
void locFstretch(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2,
                 const double& l0, const double& ks, Eigen::VectorXd& loc_f);
void locJstretch(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2,
                 const double& l0, const double& ks, Eigen::MatrixXd& loc_j);
void locFshear(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
               const double& a0, const double& ksh, Eigen::VectorXd& loc_f);
void locJshear(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
               const double& a0, const double& ksh, Eigen::MatrixXd& loc_j);
void locFbend(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
              const Eigen::Vector3d& kap0, const double& n1n20, Eigen::VectorXd& loc_f);
void locJbend(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
              const Eigen::Vector3d& kap0, const double& n1n20, Eigen::MatrixXd& loc_j);

#endif //PLATES_SHELLS_SOLVER_H
