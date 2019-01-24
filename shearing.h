#ifndef PLATES_SHELLS_SHEAR_DERIVATIVES_H
#define PLATES_SHELLS_SHEAR_DERIVATIVES_H

#include <Eigen/Dense>

/*
 *      Shearing energy
 *
 *      E_sh = 1/2 * k_sh * (phi - phi0)^2
 *
 */

class Element;

class Shearing {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Shearing(Element* ptr, double E, double nu, double area, double clen);
    void initValues();
    void locShear(Eigen::VectorXd& loc_f, Eigen::MatrixXd& loc_j);

private:
    void grad(Eigen::VectorXd& gradPhi);
    void hess(Eigen::VectorXd& gradPhi, Eigen::MatrixXd& hessPhi);

    Element* m_element;

    Eigen::Vector3d m_e1;
    Eigen::Vector3d m_e2;
    Eigen::Matrix3d m_M1;
    Eigen::Matrix3d m_M2;

    double m_phi;
    double m_ne1;
    double m_ne2;
    double m_h1;
    double m_h2;
    double m_ksh;           // Shear stiffness
};

#endif //PLATES_SHELLS_SHEAR_DERIVATIVES_H
