#include <cmath>
#include "shearing.h"
#include "element.h"
#include "node.h"

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

Shearing::Shearing(Element* ptr, double E, double nu, double area, double clen) {
    m_element = ptr;

    double G = E / (2.0 * (1.0 + nu));
    m_ksh = G * area * clen;

    initValues();
}

void Shearing::initValues() {
    m_e1 = m_element->get_node(1)->get_xyz() - m_element->get_node(2)->get_xyz();
    m_e2 = m_element->get_node(3)->get_xyz() - m_element->get_node(2)->get_xyz();

    m_ne1 = m_e1.norm();
    m_ne2 = m_e2.norm();

    m_e1 = m_e1 / m_ne1;
    m_e2 = m_e2 / m_ne2;

    double val = (m_e1 - m_e2).norm() / (m_e1 + m_e2).norm();
    m_phi = 2 * atan(val);

    m_h1 = m_ne2 * sin(m_phi);
    m_h2 = m_ne1 * sin(m_phi);

    Eigen::Matrix3d id3 = Eigen::Matrix3d::Identity(3,3);

    m_M1 = id3 - m_e1 * m_e1.transpose();
    m_M2 = id3 - m_e2 * m_e2.transpose();
}

// -----------------------------------------------------------------------

void Shearing::locShear(Eigen::VectorXd& loc_f, Eigen::MatrixXd& loc_j) {
    Eigen::VectorXd gradPhi = Eigen::VectorXd::Zero(9);
    Eigen::MatrixXd hessPhi = Eigen::MatrixXd::Zero(9, 9);

    grad(gradPhi);
    hess(gradPhi, hessPhi);

    loc_f = m_ksh * (m_phi - m_element->get_phi0()) * gradPhi;
    loc_j = m_ksh * (gradPhi * gradPhi.transpose() + (m_phi - m_element->get_phi0()) * hessPhi);
}

// -----------------------------------------------------------------------

void Shearing::grad(Eigen::VectorXd& gradPhi) {
    gradPhi.segment(0, 3) = - m_e2.transpose() * m_M1 / m_h2;
    gradPhi.segment(6, 3) = - m_e1.transpose() * m_M2 / m_h1;
    gradPhi.segment(3, 3) = - gradPhi.segment(6, 3) - gradPhi.segment(0, 3);
}

void Shearing::hess(Eigen::VectorXd& gradPhi, Eigen::MatrixXd& hessPhi) {

    double C1 = m_ne1 * cos(m_phi);
    double C2 = m_ne2 * cos(m_phi);

    Eigen::Vector3d K21 = m_e2 - m_e1 * cos(m_phi);
    Eigen::Vector3d K12 = m_e1 - m_e2 * cos(m_phi);

    Eigen::RowVector3d R1 = m_e1.transpose() * sin(m_phi);
    Eigen::RowVector3d R2 = m_e2.transpose() * sin(m_phi);

    Eigen::Matrix3d N1 = m_M1 / m_ne1;
    Eigen::Matrix3d N2 = m_M2 / m_ne2;

    Eigen::Matrix3d M11 = sin(m_phi) * m_e1 * gradPhi.segment(0, 3).transpose();
    Eigen::Matrix3d M12 = sin(m_phi) * m_e1 * gradPhi.segment(3, 3).transpose();
    Eigen::Matrix3d M13 = sin(m_phi) * m_e1 * gradPhi.segment(6, 3).transpose();
    Eigen::Matrix3d M22 = sin(m_phi) * m_e2 * gradPhi.segment(3, 3).transpose();
    Eigen::Matrix3d M33 = sin(m_phi) * m_e2 * gradPhi.segment(6, 3).transpose();


    hessPhi.block(0,0, 3,3) = - 1 / m_h2 * (M11 - N1 * cos(m_phi)) + 1 / pow(m_h2, 2) * K21 * (R1 + C1 * gradPhi.segment(0, 3).transpose());
    hessPhi.block(0,3, 3,3) = - 1 / m_h2 * (- N2 + N1 * cos(m_phi) + M12) + 1 / pow(m_h2, 2) * K21 * (- R1 + C1 * gradPhi.segment(3, 3).transpose());
    hessPhi.block(0,6, 3,3) = - 1 / m_h2 * (N2 + M13) + 1 / pow(m_h2, 2) * K21 * (C1 * gradPhi.segment(6, 3).transpose());
    hessPhi.block(3,3, 3,3) = ( 1 / m_h1 * (- N1 + N2 * cos(m_phi) + M22) + 1 / pow(m_h1, 2) * K12 * (R2 - C2 * gradPhi.segment(3, 3).transpose()) )
                              - hessPhi.block(0,3, 3,3);

    hessPhi.block(6,6, 3,3) = - 1 / m_h1 * (- N1 * cos(m_phi) + M33) + 1 / pow(m_h1, 2) * K12 * (R2 + C2 * gradPhi.segment(6, 3).transpose());
    hessPhi.block(3,6, 3,3) = - hessPhi.block(6,6, 3,3) - hessPhi.block(0,6, 3,3);

    // symmetric matrix
    hessPhi.block(3,0, 3,3) = hessPhi.block(0,3, 3,3);
    hessPhi.block(6,0, 3,3) = hessPhi.block(0,6, 3,3);
    hessPhi.block(6,3, 3,3) = hessPhi.block(3,6, 3,3);
}
