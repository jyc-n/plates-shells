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

Stretching::Stretching(Edge* ptr, double E, double D) {
    m_edge = ptr;

    m_ne0 = m_edge->get_node(2)->get_xyz() - m_edge->get_node(1)->get_xyz();
    m_len = m_ne0.norm();
    m_ne0 = m_ne0 / m_len;

    m_ks = E * M_PI/4.0 * pow(D, 2) / m_len;
}

// -----------------------------------------------------------------------
void Stretching::locStretch(Eigen::VectorXd& loc_f, Eigen::MatrixXd& loc_j) {
    Eigen::VectorXd gradLen = Eigen::VectorXd::Zero(6);
    Eigen::MatrixXd hessLen = Eigen::MatrixXd::Zero(6, 6);

    grad(gradLen);
    hess(hessLen);

    loc_f = m_ks * gradLen;
    loc_j = m_ks * hessLen;
}

// -----------------------------------------------------------------------

void Stretching::grad(Eigen::VectorXd& gradLen) {
    double deltaLen = m_len - m_edge->get_len0();
    gradLen.segment(0, 3) = - deltaLen * m_ne0;
    gradLen.segment(3, 3) = deltaLen * m_ne0;
}

void Stretching::hess(Eigen::MatrixXd& hessLen) {
    double deltaRatio = 1 - m_edge->get_len0() / m_len;
    Eigen::Matrix3d dyadicMat = m_ne0 * m_ne0.transpose();
    Eigen::Matrix3d id3 = Eigen::Matrix3d::Identity(3,3);

    hessLen.block(0,0, 3,3) = dyadicMat + deltaRatio * (id3 - dyadicMat.transpose());
    hessLen.block(0,3, 3,3) = - hessLen.block(0,0, 3,3);
    hessLen.block(3,3, 3,3) = hessLen.block(0,0, 3,3);

    // symmetric matrix
    hessLen.block(3,0, 3,3) = hessLen.block(0,3, 3,3);
}