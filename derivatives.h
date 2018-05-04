#ifndef PLATES_SHELLS_DERIVATIVES_H
#define PLATES_SHELLS_DERIVATIVES_H

#include <Eigen/Dense>

/*
 *      Stretch energy
 *
 *      E_s = k_s * l0^2 * (l/l0 - 1)^2
 *
 *      l0: original length, regarded as constant
 *      l: length after deformation
 *
 *      l = sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)
 */

// 1st derivatives of stretch energy
double dEs_dx1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double dEs_dx2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double dEs_dy1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double dEs_dy2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double dEs_dz1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double dEs_dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);

// 2nd derivatives of stretch energy
double ddEs_dx1dx1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double ddEs_dx1dx2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double ddEs_dx1dy1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double ddEs_dx1dy2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double ddEs_dx1dz1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double ddEs_dx1dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);

double ddEs_dx2dx2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double ddEs_dx2dy1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double ddEs_dx2dy2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double ddEs_dx2dz1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double ddEs_dx2dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);

double ddEs_dy1dy1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double ddEs_dy1dy2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double ddEs_dy1dz1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double ddEs_dy1dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);

double ddEs_dy2dy2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double ddEs_dy2dz1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double ddEs_dy2dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);

double ddEs_dz1dz1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double ddEs_dz1dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);
double ddEs_dz2dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks);

/*
*      Shear energy
*
*      E_sh = 1/2 * k_sh * (psi - psi0)^2
*
*/

/*
 *      Bending energy
 *
 *      E_b = 1/2 * k_b * (theta - theta0)^2
 *
 */

#endif //PLATES_SHELLS_DERIVATIVES_H
