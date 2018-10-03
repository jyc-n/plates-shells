#ifndef PLATES_SHELLS_BENDING_DERIVATIVES_H
#define PLATES_SHELLS_BENDING_DERIVATIVES_H

#include <Eigen/Dense>

/*
 *      Bending energy
 *
 *      E_b = k_b * ( 3 * l0^2 )/A * ( phi(theta) - phi0 )^2
 *
 *      grad E = zeta(theta, theta0) * grad theta
 *      Hess E = zeta(theta, theta0) * Hess theta + xi(theta, theta0) * grad theta * grad theta ^T
 *
 */

class Hinge;

class Bending {
public:
    Bending(Hinge* ptr);
    void init();
    void locBend(Eigen::VectorXd& loc_f, Eigen::MatrixXd& loc_j);

    double m_coeff;

private:
    void perturb(int pos, double val);
    void recover(int pos);
    double getAngle();
    double getEnergy();

    void grad();
    void hess();

    bool READY;

    Hinge* m_hinge;

    Eigen::VectorXd m_q;
    Eigen::VectorXd m_qCurrent;
    Eigen::VectorXd m_gradE;
    Eigen::MatrixXd m_hessE;

    double m_phi0;
    double m_phiPerturbed;

    const double delta;
};


#endif //PLATES_SHELLS_BENDING_DERIVATIVES_H
