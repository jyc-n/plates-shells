#ifndef PLATES_SHELLS_BENDING_DERIVATIVES_H
#define PLATES_SHELLS_BENDING_DERIVATIVES_H

#include <Eigen/Dense>

class Hinge;

/*
 *      Bending energy
 *
 *      E_b = k_b * ( 3 * l0^2 )/A * ( phi(theta) - phi0 )^2
 *
 *      - grad E = zeta(theta, theta0) * grad theta
 *      - Hess E = zeta(theta, theta0) * Hess theta + xi(theta, theta0) * grad theta ^T * grad theta
 *
 */

class Bending {
public:
    Bending(Hinge* ptr);

    void initValues();
    void locBend(Eigen::VectorXd& loc_f, Eigen::MatrixXd& loc_j);

private:
    void psi();
    void zeta();
    void xi();
    Eigen::Matrix3d s(Eigen::Matrix3d& mat);

    void grad(Eigen::VectorXd& gradTheta);
    void hess(Eigen::MatrixXd& hessTheta);


    Hinge* m_hinge;

    Eigen::Vector3d m_e0;
    Eigen::Vector3d m_e1;
    Eigen::Vector3d m_e2;
    Eigen::Vector3d m_e3;
    Eigen::Vector3d m_e4;

    Eigen::Vector3d m_nn1;
    Eigen::Vector3d m_nn2;
    Eigen::Vector3d m_m1;
    Eigen::Vector3d m_m2;
    Eigen::Vector3d m_m3;
    Eigen::Vector3d m_m4;
    Eigen::Vector3d m_m01;
    Eigen::Vector3d m_m02;

    double m_theta;
    double m_alpha1;
    double m_alpha2;
    double m_alpha3;
    double m_alpha4;

    double m_psi;
    double m_zeta;
    double m_xi;

    double m_h1;
    double m_h2;
    double m_h3;
    double m_h4;
    double m_h01;
    double m_h02;
};


#endif //PLATES_SHELLS_BENDING_DERIVATIVES_H
