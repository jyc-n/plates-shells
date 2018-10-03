#include <cmath>
#include "shearing.h"
#include "element.h"
#include "node.h"
#include <iostream>

/*
 *                x3
 *                 ^
 *               / |
 *             /   |
 *           /     |    e2
 *         /       |
 *       /         |
 *   x1 <----------  x2
 *          e1
 *
 */

// -----------------------------------------------------------------------

Shearing::Shearing(Element* ptr)
  : delta(1e-6)
{
    m_element = ptr;
    READY = false;
}

void Shearing::init() {
    m_phi0 = m_element->get_phi0();

    // dof vector x1 y1 z1 x2 y2 z2 x3 y3 z3
    m_q = Eigen::VectorXd::Zero(9);

    m_q.segment(0,3) = m_element->get_node(1)->get_xyz();
    m_q.segment(3,3) = m_element->get_node(2)->get_xyz();
    m_q.segment(6,3) = m_element->get_node(3)->get_xyz();

    m_qCurrent = m_q;

    m_gradE = Eigen::VectorXd::Zero(9);
    m_hessE = Eigen::MatrixXd::Zero(9, 9);

    READY = true;
}

// -----------------------------------------------------------------------

void Shearing::perturb(int pos, double val) {
    m_qCurrent(pos) += val;
    m_phiPerturbed = getAngle();
}

void Shearing::recover(int pos) {
    m_qCurrent(pos) = m_q(pos);
    m_phiPerturbed = getAngle();
}

double Shearing::getAngle() {
    Eigen::Vector3d e1 = Eigen::Vector3d::Zero();
    Eigen::Vector3d e2 = Eigen::Vector3d::Zero();

    e1 = m_qCurrent.segment(0,3) - m_qCurrent.segment(3,3);
    e2 = m_qCurrent.segment(6,3) - m_qCurrent.segment(3,3);

    double ne1 = e1.norm();
    double ne2 = e2.norm();

    e1 = e1 / ne1;
    e2 = e2 / ne2;

    double val = (e1 - e2).norm() / (e1 + e2).norm();
    double phi = 2 * atan(val);

    return phi;
}

double Shearing::getEnergy() {
    return m_coeff * pow(m_phiPerturbed - m_phi0, 2);
}

// -----------------------------------------------------------------------

void Shearing::locShear(Eigen::VectorXd& loc_f, Eigen::MatrixXd& loc_j) {
    if (!READY)
        throw "must initialize before solving";

    grad();
    hess();

    loc_f = m_gradE;
    loc_j = m_hessE;
}

void Shearing::grad() {
    double EshHi = 0, EshLo = 0;
    for (int i = 0; i < 9; i++) {
        perturb(i, delta);
        EshHi = getEnergy();
        recover(i);

        perturb(i, -delta);
        EshLo = getEnergy();
        recover(i);

        m_gradE(i) = (EshHi - EshLo) / (2 * delta);
    }
}

void Shearing::hess() {
    double Esh_i_Hi_j_Hi = 0, Esh_i_Hi_j_Lo = 0, Esh_i_Lo_j_Hi = 0, Esh_i_Lo_j_Lo = 0;
    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 9; j++) {
            perturb(i, delta);
            perturb(j, delta);
            Esh_i_Hi_j_Hi = getEnergy();
            recover(i);
            recover(j);

            perturb(i, delta);
            perturb(j, -delta);
            Esh_i_Hi_j_Lo = getEnergy();
            recover(i);
            recover(j);

            perturb(i, -delta);
            perturb(j, delta);
            Esh_i_Lo_j_Hi = getEnergy();
            recover(i);
            recover(j);

            perturb(i, -delta);
            perturb(j, -delta);
            Esh_i_Lo_j_Lo = getEnergy();
            recover(i);
            recover(j);

            m_hessE(i, j) = (Esh_i_Hi_j_Hi - Esh_i_Hi_j_Lo - Esh_i_Lo_j_Hi + Esh_i_Lo_j_Lo) / (4 * pow(delta, 2));
        }
    }
}
