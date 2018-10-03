#ifndef PLATES_SHELLS_DERIVATIVES_H
#define PLATES_SHELLS_DERIVATIVES_H

#include <Eigen/Dense>

/*
 *      Stretch energy (numerical version)
 *
 *      E_s = 1/2 * k_s * l0^2 * (l/l0 - 1)^2
 *
 *      l0: original length, regarded as constant
 *      l: length after deformation
 *
 */

class Edge;

class Stretching {
public:
    Stretching(Edge* ptr);
    void init();
    void locStretch(Eigen::VectorXd& loc_f, Eigen::MatrixXd& loc_j);

    double m_ks;

private:
    void perturb(int pos, double val);
    void recover(int pos);
    double getLen();
    double getEnergy();

    void grad();
    void hess();

    bool READY;
    double m_coeff;

    Edge* m_edge;

    Eigen::VectorXd m_q;
    Eigen::VectorXd m_qCurrent;
    Eigen::VectorXd m_gradE;
    Eigen::MatrixXd m_hessE;

    double m_l0;
    double m_lPerturbed;

    const double delta;
};

#endif //PLATES_SHELLS_DERIVATIVES_H
