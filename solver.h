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

    // main solver function
    void Solve();

    // subroutine
    void updateResidual(Eigen::VectorXd& qn, Eigen::VectorXd& qnew, Eigen::VectorXd& vel, Eigen::VectorXd& res_f);
    void updateJacobian(Eigen::MatrixXd& mat_j);
    Eigen::VectorXd calcVel(double dt, const Eigen::VectorXd& qcurr, const Eigen::VectorXd& qnew);
    Eigen::VectorXd calcDofnew(Eigen::VectorXd& qn, const Eigen::VectorXd& temp_f, const Eigen::MatrixXd& temp_j);

    // helper functions
    Eigen::VectorXd unconsVec(const Eigen::VectorXd& vec);
    Eigen::MatrixXd unconsMat(const Eigen::MatrixXd& mat);
    void updateNodes(Eigen::VectorXd& qnew);
    void resetMarks();

    void calculate_fstretch(Element& el, Eigen::VectorXd& dEdq);
    void calculate_fshear(Element& el, Eigen::VectorXd& dEdq);
    void calculate_fbend(Element& el, Eigen::VectorXd& dEdq);
    void calculate_jstretch(Element& el, Eigen::MatrixXd& ddEddq);
    void calculate_jshear(Element& el, Eigen::MatrixXd& ddEddq);
    void calculate_jbend(Element& el, Eigen::MatrixXd& ddEddq);


private:
    Parameters* m_SimPar;
    Geometry*   m_SimGeo;
    Boundary*   m_SimBC;
};

#endif //PLATES_SHELLS_SOLVER_H
