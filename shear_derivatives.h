#ifndef PLATES_SHELLS_SHEAR_DERIVATIVES_H
#define PLATES_SHELLS_SHEAR_DERIVATIVES_H

#include <Eigen/Dense>

/*
 *      Shear energy
 *
 *      E_sh = 1/2 * k_sh * A0 * (A/A0 - 1)^2
 *
 *      A0: original area
 *      A: current area
 *
 *      A = norm(e1 x e2)
 *      A0 = norm(e1) * norm(e2)
 */

// 1st derivatives of shear energy
double dEsh_dx1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                const double& a0, const double& ksh);
double dEsh_dx2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                const double& a0, const double& ksh);
double dEsh_dx3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                const double& a0, const double& ksh);
double dEsh_dy1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                const double& a0, const double& ksh);
double dEsh_dy2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                const double& a0, const double& ksh);
double dEsh_dy3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                const double& a0, const double& ksh);
double dEsh_dz1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                const double& a0, const double& ksh);
double dEsh_dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                const double& a0, const double& ksh);
double dEsh_dz3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                const double& a0, const double& ksh);

// 2nd derivatives of shear energy
double ddEsh_dx1dx1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx1dx2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx1dx3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx1dy1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx1dy2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx1dy3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx1dz1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx1dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx1dz3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);

double ddEsh_dx2dx2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx2dx3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx2dy1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx2dy2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx2dy3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx2dz1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx2dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx2dz3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);

double ddEsh_dx3dx3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx3dy1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx3dy2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx3dy3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx3dz1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx3dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dx3dz3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);

double ddEsh_dy1dy1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dy1dy2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dy1dy3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dy1dz1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dy1dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dy1dz3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);

double ddEsh_dy2dy2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dy2dy3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dy2dz1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dy2dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dy2dz3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);

double ddEsh_dy3dy3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dy3dz1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dy3dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dy3dz3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);

double ddEsh_dz1dz1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dz1dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dz1dz3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);

double ddEsh_dz2dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dz2dz3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);
double ddEsh_dz3dz3(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
                    const double& a0, const double& ksh);

#endif //PLATES_SHELLS_SHEAR_DERIVATIVES_H
