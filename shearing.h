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
    Shearing(Element* ptr);
    void init();
    void locShear(Eigen::VectorXd& loc_f, Eigen::MatrixXd& loc_j);

    double m_coeff;

private:
    void perturb(int pos, double val);
    void recover(int pos);
    double getAngle();
    double getEnergy();

    void grad();
    void hess();

    bool READY;

    Element* m_element;

    Eigen::VectorXd m_q;
    Eigen::VectorXd m_qCurrent;
    Eigen::VectorXd m_gradE;
    Eigen::MatrixXd m_hessE;

    double m_phi0;
    double m_phiPerturbed;

    const double delta;
};

#endif //PLATES_SHELLS_SHEAR_DERIVATIVES_H
