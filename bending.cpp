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

Bending::Bending(Hinge* ptr) {
    m_hinge = ptr;

    m_e0 = m_hinge->get_node(1)->get_xyz() - m_hinge->get_node(0)->get_xyz();
    m_e1 = m_hinge->get_node(2)->get_xyz() - m_hinge->get_node(0)->get_xyz();
    m_e2 = m_hinge->get_node(3)->get_xyz() - m_hinge->get_node(0)->get_xyz();
    m_e3 = m_hinge->get_node(2)->get_xyz() - m_hinge->get_node(1)->get_xyz();
    m_e4 = m_hinge->get_node(3)->get_xyz() - m_hinge->get_node(1)->get_xyz();
    initValues();
}

void Bending::initValues() {
    m_alpha1 = abs(acos(m_e0.dot(m_e1) / (m_e0.norm() * m_e1.norm())));
    m_alpha2 = abs(acos(m_e0.dot(m_e2) / (m_e0.norm() * m_e2.norm())));
    m_alpha3 = M_PI - abs(acos(m_e0.dot(m_e3) / (m_e0.norm() * m_e3.norm())));
    m_alpha4 = M_PI - abs(acos(m_e0.dot(m_e4) / (m_e0.norm() * m_e4.norm())));

    m_nn1 = m_e0.cross(m_e3);
    m_nn1 = m_nn1 / (m_nn1.norm());
    m_nn2 = -m_e0.cross(m_e4);
    m_nn2 = m_nn2 / (m_nn2.norm());

    m_m1 = (m_nn1).cross(m_e1 / m_e1.norm());
    m_m2 = (m_e2 / m_e2.norm()).cross(m_nn2);
    m_m3 = (m_e3 / m_e3.norm()).cross(m_nn1);
    m_m4 = (m_e4 / m_e4.norm()).cross(m_nn2);
    m_m01 = (m_e0 / m_e0.norm()).cross(m_nn1);
    m_m02 = (m_nn2).cross(m_e0 / m_e0.norm());

    double val = (m_nn1 - m_nn2).norm() / (m_nn1 + m_nn2).norm();
    m_theta = 2 * atan(val);

    m_h1 = m_e0.norm() * sin(m_alpha1);
    m_h2 = m_e0.norm() * sin(m_alpha2);
    m_h3 = m_e0.norm() * sin(m_alpha3);
    m_h4 = m_e0.norm() * sin(m_alpha4);
    m_h01 = m_e1.norm() * sin(m_alpha1);
    m_h02 = m_e2.norm() * sin(m_alpha2);

    psi();
    zeta();
    xi();
}

// -----------------------------------------------------------------------

void Bending::locBend(Eigen::VectorXd &loc_f, Eigen::MatrixXd &loc_j) {
    Eigen::VectorXd gradTheta = Eigen::VectorXd::Zero(12);
    Eigen::MatrixXd hessTheta = Eigen::MatrixXd::Zero(12, 12);

    grad(gradTheta);
    hess(hessTheta);

    loc_f = m_zeta * gradTheta;
    loc_j = m_zeta * hessTheta + m_xi * gradTheta * gradTheta.transpose();
}

// -----------------------------------------------------------------------

void Bending::psi() {
    m_psi = tan(m_theta / 2.0);
}

void Bending::zeta() {
    m_zeta = (2.0 * m_hinge->m_k * (m_psi - m_hinge->m_psi0) * (1 + pow(m_psi, 2)));
}

void Bending::xi() {
    m_xi = (m_hinge->m_k * (1 + pow(m_psi, 2)) * (2 * (m_psi - m_hinge->m_psi0) * m_psi + (1 + pow(m_psi, 2))));
}

Eigen::Matrix3d Bending::s(Eigen::Matrix3d &mat) {
    return (mat + mat.transpose());
}

void Bending::grad(Eigen::VectorXd& gradTheta) {
    gradTheta.segment(0, 3) = cos(m_alpha3) * m_nn1 / m_h3 + cos(m_alpha4) * m_nn2 / m_h4;
    gradTheta.segment(3, 3) = cos(m_alpha1) * m_nn1 / m_h1 + cos(m_alpha2) * m_nn2 / m_h2;
    gradTheta.segment(6, 3) = - m_nn1 / m_h01;
    gradTheta.segment(9, 3) = - m_nn2 / m_h02;
}

void Bending::hess(Eigen::MatrixXd& hessTheta) {
    Eigen::Matrix3d M331 = cos(m_alpha3) / (m_h3 * m_h3) * m_m3 * m_nn1.transpose();
    Eigen::Matrix3d M311 = cos(m_alpha3) / (m_h3 * m_h1) * m_m1 * m_nn1.transpose();
    Eigen::Matrix3d M131 = cos(m_alpha1) / (m_h1 * m_h3) * m_m3 * m_nn1.transpose();
    Eigen::Matrix3d M3011 = cos(m_alpha3) / (m_h3 * m_h01) * m_m01 * m_nn1.transpose();
    Eigen::Matrix3d M111 = cos(m_alpha1) / (m_h1 * m_h1) * m_m1 * m_nn1.transpose();
    Eigen::Matrix3d M1011 = cos(m_alpha1) / (m_h1 * m_h01) * m_m01 * m_nn1.transpose();

    Eigen::Matrix3d M442 = cos(m_alpha4) / (m_h4 * m_h4) * m_m4 * m_nn2.transpose();
    Eigen::Matrix3d M422 = cos(m_alpha4) / (m_h4 * m_h2) * m_m2 * m_nn2.transpose();
    Eigen::Matrix3d M242 = cos(m_alpha2) / (m_h2 * m_h4) * m_m4 * m_nn2.transpose();
    Eigen::Matrix3d M4022 = cos(m_alpha4) / (m_h4 * m_h02) * m_m02 * m_nn2.transpose();
    Eigen::Matrix3d M222 = cos(m_alpha2) / (m_h2 * m_h2) * m_m2 * m_nn2.transpose();
    Eigen::Matrix3d M2022 = cos(m_alpha2) / (m_h2 * m_h02) * m_m02 * m_nn2.transpose();

    Eigen::Matrix3d B1 = 1 / (pow(m_e0.norm(), 2)) * m_nn1 * m_m01.transpose();
    Eigen::Matrix3d B2 = 1 / (pow(m_e0.norm(), 2)) * m_nn2 * m_m02.transpose();

    Eigen::Matrix3d N13 = 1 / (m_h01 * m_h3) * m_nn1 * m_m3.transpose();
    Eigen::Matrix3d N24 = 1 / (m_h02 * m_h4) * m_nn2 * m_m4.transpose();
    Eigen::Matrix3d N11 = 1 / (m_h01 * m_h1) * m_nn1 * m_m1.transpose();
    Eigen::Matrix3d N22 = 1 / (m_h02 * m_h2) * m_nn2 * m_m2.transpose();
    Eigen::Matrix3d N101 = 1 / (m_h01 * m_h01) * m_nn1 * m_m01.transpose();
    Eigen::Matrix3d N202 = 1 / (m_h02 * m_h02) * m_nn2 * m_m02.transpose();

    hessTheta.block(0,0, 3,3) = s(M331) - B1 + s(M442) - B2;
    hessTheta.block(0,3, 3,3) = M331 + M131.transpose() + B1 + M422 + M242.transpose() + B2;
    hessTheta.block(0,6, 3,3) = M3011 - N13;
    hessTheta.block(0,9, 3,3) = M4022 - N24;
    hessTheta.block(3,3, 3,3) = s(M111) - B1 + s(M222) - B2;
    hessTheta.block(3,6, 3,3) = M1011 - N11;
    hessTheta.block(3,9, 3,3) = M2022 - N22;
    hessTheta.block(6,6, 3,3) = -s(N101);
    hessTheta.block(9,9, 3,3) = -s(N202);

    // symmetric matrix
    hessTheta.block(3,0, 3,3) = hessTheta.block(0,3, 3,3);
    hessTheta.block(6,0, 3,3) = hessTheta.block(0,6, 3,3);
    hessTheta.block(9,0, 3,3) = hessTheta.block(0,9, 3,3);
    hessTheta.block(6,3, 3,3) = hessTheta.block(3,6, 3,3);
    hessTheta.block(9,3, 3,3) = hessTheta.block(9,6, 3,3);
}