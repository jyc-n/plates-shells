#ifndef PLATES_SHELLS_DERIVATIVES_H
#define PLATES_SHELLS_DERIVATIVES_H

#include <Eigen/Dense>

/*
 *      Stretch energy
 *
 *      E_s = k_s * l0^2 * (l/l0 - 1)^2
 *
 *      l0: original length, regarded as constant
 *      l: length after deformation
 *
 */

class Edge;

class Stretching {
public:
    Stretching(Edge* ptr);

    void locStretch(Eigen::VectorXd& loc_f, Eigen::MatrixXd& loc_j);

private:
    void grad(Eigen::VectorXd& gradLen);
    void hess(Eigen::MatrixXd& hessLen);

    Edge* m_edge;

    Eigen::Vector3d m_ne0;
    double m_len;
};

#endif //PLATES_SHELLS_DERIVATIVES_H
