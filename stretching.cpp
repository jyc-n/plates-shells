#include <cmath>
#include "stretching.h"
#include "edge.h"
#include "node.h"

/*
 *     x1                x2
 *     ------------------>
 *              e0
 */

// -----------------------------------------------------------------------
Stretching::Stretching(Edge* ptr)
  : delta(1e-6)
{
    m_edge = ptr;
    READY = false;
}

void Stretching::init() {
    m_l0 = m_edge->get_len0();
    m_coeff = 0.5 * m_ks * pow(m_l0, 2);

    // dof vector x1 y1 z1 x2 y2 z2
    m_q = Eigen::VectorXd::Zero(6);

    m_q.segment(0,3) = m_edge->get_node(1)->get_xyz();
    m_q.segment(3,3) = m_edge->get_node(2)->get_xyz();

    m_qCurrent = m_q;

    m_gradE = Eigen::VectorXd::Zero(6);
    m_hessE = Eigen::MatrixXd::Zero(6, 6);

    READY = true;
}

// -----------------------------------------------------------------------
void Stretching::perturb(int pos, double val) {
    m_qCurrent(pos) += val;
    m_lPerturbed = getLen();
}

void Stretching::recover(int pos) {
    m_qCurrent(pos) = m_q(pos);
    m_lPerturbed = getLen();
}

double Stretching::getLen() {
    return sqrt(pow(m_qCurrent(0) - m_qCurrent(3), 2)
              + pow(m_qCurrent(1) - m_qCurrent(4), 2)
              + pow(m_qCurrent(2) - m_qCurrent(5), 2));
}

double Stretching::getEnergy() {
    return m_coeff * pow(m_lPerturbed/m_l0 - 1, 2);
}

// -----------------------------------------------------------------------
void Stretching::locStretch(Eigen::VectorXd& loc_f, Eigen::MatrixXd& loc_j) {
    if (!READY)
        throw "must initialize before solving";

    grad();
    hess();

    loc_f = m_gradE;
    loc_j = m_hessE;
}

void Stretching::grad() {
    double EsHi = 0, EsLo = 0;
    for (int i = 0; i < 6; i++) {
        perturb(i, delta);
        EsHi = getEnergy();
        recover(i);

        perturb(i, -delta);
        EsLo = getEnergy();
        recover(i);

        m_gradE(i) = (EsHi - EsLo) / (2 * delta);
    }
}

void Stretching::hess() {
    double Es_i_Hi_j_Hi = 0, Es_i_Hi_j_Lo = 0, Es_i_Lo_j_Hi = 0, Es_i_Lo_j_Lo = 0;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            perturb(i, delta);
            perturb(j, delta);
            Es_i_Hi_j_Hi = getEnergy();
            recover(i);
            recover(j);

            perturb(i, delta);
            perturb(j, -delta);
            Es_i_Hi_j_Lo = getEnergy();
            recover(i);
            recover(j);

            perturb(i, -delta);
            perturb(j, delta);
            Es_i_Lo_j_Hi = getEnergy();
            recover(i);
            recover(j);

            perturb(i, -delta);
            perturb(j, -delta);
            Es_i_Lo_j_Lo = getEnergy();
            recover(i);
            recover(j);

            m_hessE(i, j) = (Es_i_Hi_j_Hi - Es_i_Hi_j_Lo - Es_i_Lo_j_Hi + Es_i_Lo_j_Lo) / (4 * pow(delta, 2));
        }
    }
}