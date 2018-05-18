#ifndef PLATES_SHELLS_BENDING_DERIVATIVES_H
#define PLATES_SHELLS_BENDING_DERIVATIVES_H

#include <Eigen/Dense>

/*
 *      Bending energy
 *
 *      E_b = k_b * ( 3 * l0^2 )/A * ( phi(theta) - phi0 )^2
 *
 *      l0: original hinge length
 *      A: sum of area
 *      theta: current dihedral angle
 *      theta0: original dihedral angle
 *
 *      l0 = norm(e0)                                           constant
 *      A = A1 + A2                                             constant
 *      phi0 = 2 * tan( theta0/2 )                              constant
 *      phi(theta) = 2 * tan( theta/2 )
 *      tan(theta/2) = norm(nn1 - nn2) / norm(nn1 + nn2)
 *      nn1 = n1 / norm(n1)
 *      nn2 = n2 / norm(n2)
 *      n1 = e0 x e1
 *      n2 = e2 x e0
 *
 */

// 1st derivatives of bending energy
double dEb_dx0(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
               const double& l0, const double& kb, const double& a0, const double& phi0);
double dEb_dx1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
               const double& l0, const double& kb, const double& a0, const double& phi0);
double dEb_dx2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
               const double& l0, const double& kb, const double& a0, const double& phi0);
double dEb_dx3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
               const double& l0, const double& kb, const double& a0, const double& phi0);
double dEb_dy0(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
               const double& l0, const double& kb, const double& a0, const double& phi0);
double dEb_dy1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
               const double& l0, const double& kb, const double& a0, const double& phi0);
double dEb_dy2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
               const double& l0, const double& kb, const double& a0, const double& phi0);
double dEb_dy3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
               const double& l0, const double& kb, const double& a0, const double& phi0);
double dEb_dz0(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
               const double& l0, const double& kb, const double& a0, const double& phi0);
double dEb_dz1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
               const double& l0, const double& kb, const double& a0, const double& phi0);
double dEb_dz2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
               const double& l0, const double& kb, const double& a0, const double& phi0);
double dEb_dz3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
               const double& l0, const double& kb, const double& a0, const double& phi0);

// 2nd derivatives of bending energy
double ddEb_dx0dx0(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx0dx1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx0dx2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx0dx3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx0dy0(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx0dy1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx0dy2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx0dy3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx0dz0(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx0dz1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx0dz2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx0dz3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);

double ddEb_dx1dx1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx1dx2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx1dx3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx1dy0(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx1dy1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx1dy2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx1dy3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx1dz0(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx1dz1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx1dz2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx1dz3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);

double ddEb_dx2dx2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx2dx3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx2dy0(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx2dy1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx2dy2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx2dy3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx2dz0(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx2dz1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx2dz2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx2dz3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);

double ddEb_dx3dx3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx3dy0(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx3dy1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx3dy2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx3dy3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx3dz0(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx3dz1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx3dz2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dx3dz3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);

double ddEb_dy0dy0(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy0dy1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy0dy2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy0dy3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy0dz0(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy0dz1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy0dz2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy0dz3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);

double ddEb_dy1dy1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy1dy2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy1dy3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy1dz0(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy1dz1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy1dz2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy1dz3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);

double ddEb_dy2dy2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy2dy3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy2dz0(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy2dz1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy2dz2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy2dz3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);

double ddEb_dy3dy3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy3dz0(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy3dz1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy3dz2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dy3dz3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);

double ddEb_dz0dz0(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dz0dz1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dz0dz2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dz0dz3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);

double ddEb_dz1dz1(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dz1dz2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dz1dz3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);

double ddEb_dz2dz2(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dz2dz3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);
double ddEb_dz3dz3(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                   const double& l0, const double& kb, const double& a0, const double& phi0);

#endif //PLATES_SHELLS_BENDING_DERIVATIVES_H
