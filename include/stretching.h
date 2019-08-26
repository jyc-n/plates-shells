#ifndef PLATES_SHELLS_DERIVATIVES_H
#define PLATES_SHELLS_DERIVATIVES_H

#include <Eigen/Dense>

/*
 *      Stretch energy
 *
 *      E_s = 1/2 * k_s * (l - l0)^2
 *
 *      l0: original length, regarded as constant
 *      l: length after deformation
 *
 */

class Edge;

class Stretching {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Stretching(Edge* ptr, double E, double T);

    void locStretch(Eigen::VectorXd& loc_f, Eigen::MatrixXd& loc_j);

private:
    void grad(Eigen::VectorXd& gradLen);
    void hess(Eigen::MatrixXd& hessLen);

    Edge* m_edge;

    Eigen::Vector3d m_ne0;
    double m_len;           // Current length of edge
    double m_ks;            // Stretching stiffness
};

#endif //PLATES_SHELLS_DERIVATIVES_H
