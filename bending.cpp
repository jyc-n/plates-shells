#include <iostream>
#include <cmath>
#include "bending.h"
#include "hinge.h"
#include "node.h"

//         x2
//         /\
//        /  \
//     e1/    \e3
//      /  t0  \
//     /        \
//    /    e0    \
//  x0------------x1
//    \          /
//     \   t1   /
//      \      /
//     e2\    /e4
//        \  /
//         \/
//         x3
//
// Edge orientation: e0,e1,e2 point away from x0
//                      e3,e4 point away from x1

// -----------------------------------------------------------------------

Bending::Bending(Hinge* ptr)
  : delta(1e-6)
{
    m_hinge = ptr;
    READY = false;
}

void Bending::init() {
    m_phi0 = 2 * m_hinge->m_psi0;

    // dof vector x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4
    m_q = Eigen::VectorXd::Zero(12);

    m_q.segment(0,3) = m_hinge->get_node(0)->get_xyz();
    m_q.segment(3,3) = m_hinge->get_node(1)->get_xyz();
    m_q.segment(6,3) = m_hinge->get_node(2)->get_xyz();
    m_q.segment(9,3) = m_hinge->get_node(3)->get_xyz();

    m_qCurrent = m_q;

    m_gradE = Eigen::VectorXd::Zero(12);
    m_hessE = Eigen::MatrixXd::Zero(12, 12);

    READY = true;
}

// -----------------------------------------------------------------------

void Bending::perturb(int pos, double val) {
    m_qCurrent(pos) += val;
}

void Bending::recover(int pos) {
    m_qCurrent(pos) = m_q(pos);
}

double Bending::getAngle() {
    Eigen::Vector3d e0 = Eigen::Vector3d::Zero();
    Eigen::Vector3d e3 = Eigen::Vector3d::Zero();
    Eigen::Vector3d e4 = Eigen::Vector3d::Zero();
    Eigen::Vector3d nn1 = Eigen::Vector3d::Zero();
    Eigen::Vector3d nn2 = Eigen::Vector3d::Zero();

    e0 = m_qCurrent.segment(3,3) - m_qCurrent.segment(0,3);
    e3 = m_qCurrent.segment(6,3) - m_qCurrent.segment(3,3);
    e4 = m_qCurrent.segment(9,3) - m_qCurrent.segment(3,3);

    nn1 = e0.cross(e3);
    nn1 = nn1 / (nn1.norm());
    nn2 = -e0.cross(e4);
    nn2 = nn2 / (nn2.norm());

    double val = (nn1 - nn2).norm() / (nn1 + nn2).norm();
    double phi = 2 * atan(val);

    return phi;
}

double Bending::getEnergy() {
    m_phiPerturbed = getAngle();
    return m_coeff * pow(m_phiPerturbed - m_phi0, 2);
}

// -----------------------------------------------------------------------

void Bending::locBend(Eigen::VectorXd &loc_f, Eigen::MatrixXd &loc_j){
    if (!READY)
        throw "must initialize before solving";

    grad();
    hess();

    loc_f = m_gradE;
    loc_j = m_hessE;
}

void Bending::grad() {
    double EshHi = 0, EshLo = 0;
    for (int i = 0; i < 12; i++) {
        perturb(i, delta);
        EshHi = getEnergy();
        recover(i);

        perturb(i, -delta);
        EshLo = getEnergy();
        recover(i);

        m_gradE(i) = (EshHi - EshLo) / (2 * delta);
    }
}

void Bending::hess() {
    double Esh_i_Hi_j_Hi = 0, Esh_i_Hi_j_Lo = 0, Esh_i_Lo_j_Hi = 0, Esh_i_Lo_j_Lo = 0;
    for (int i = 0; i < 12; i++) {
        for (int j = 0; j < 12; j++) {
            if (i > j) {
                m_hessE(i, j) = m_hessE(j, i);
                continue;
            }

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