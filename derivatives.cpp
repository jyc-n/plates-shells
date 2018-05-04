<<<<<<< HEAD
//
// Created by RLA on 2018/5/2.
//

#include "derivatives.h"
=======
#include <cmath>
#include "derivatives.h"

/*
 *  ==========================================
 *      Stretch Energy, 1st derivatives
 *  ==========================================
 */

double dEs_dx1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (2*ks*l0*(v1(0) - v2(0))*(-1 + sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2))/l0))/
           sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2));
}
double dEs_dx2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (-2*ks*l0*(v1(0) - v2(0))*(-1 + sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2))/l0))/
           sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2));
}
double dEs_dy1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (2*ks*l0*(v1(1) - v2(1))*(-1 + sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2))/l0))/
           sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2));
}
double dEs_dy2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (-2*ks*l0*(v1(1) - v2(1))*(-1 + sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2))/l0))/
           sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2));
}
double dEs_dz1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (2*ks*l0*(-1 + sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2))/l0)* (v1(2) - v2(2)))/
           sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2));
}
double dEs_dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (-2*ks*l0*(-1 + sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2))/l0)* (v1(2) - v2(2)))/
           sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2));
}

/*
 *  ==========================================
 *      Stretch Energy, 2nd derivatives
 *  ==========================================
 */

double ddEs_dx1dx1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (-2*ks*l0*pow(v1(0) - v2(0),2)*
            (-1 + sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                       pow(v1(2) - v2(2),2))/l0))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) +
           (2*ks*pow(v1(0) - v2(0),2))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2)) +
           (2*ks*l0*(-1 + sqrt(pow(v1(0) - v2(0),2) +
                               pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2))/l0))/
           sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                pow(v1(2) - v2(2),2));
}
double ddEs_dx1dx2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (2*ks*l0*pow(v1(0) - v2(0),2)*
            (-1 + sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                       pow(v1(2) - v2(2),2))/l0))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) -
           (2*ks*pow(v1(0) - v2(0),2))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2)) -
           (2*ks*l0*(-1 + sqrt(pow(v1(0) - v2(0),2) +
                               pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2))/l0))/
           sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                pow(v1(2) - v2(2),2));
}
double ddEs_dx1dy1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (-2*ks*l0*(v1(0) - v2(0))*(v1(1) - v2(1))*
            (-1 + sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                       pow(v1(2) - v2(2),2))/l0))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) +
           (2*ks*(v1(0) - v2(0))*(v1(1) - v2(1)))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2));
}
double ddEs_dx1dy2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (2*ks*l0*(v1(0) - v2(0))*(v1(1) - v2(1))*
            (-1 + sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                       pow(v1(2) - v2(2),2))/l0))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) -
           (2*ks*(v1(0) - v2(0))*(v1(1) - v2(1)))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2));
}
double ddEs_dx1dz1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (-2*ks*l0*(v1(0) - v2(0))*(-1 +
                                      sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                                           pow(v1(2) - v2(2),2))/l0)*(v1(2) - v2(2)))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) +
           (2*ks*(v1(0) - v2(0))*(v1(2) - v2(2)))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2));
}
double ddEs_dx1dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (2*ks*l0*(v1(0) - v2(0))*(-1 +
                                     sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                                          pow(v1(2) - v2(2),2))/l0)*(v1(2) - v2(2)))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) -
           (2*ks*(v1(0) - v2(0))*(v1(2) - v2(2)))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2));
}

double ddEs_dx2dx2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (-2*ks*l0*pow(v1(0) - v2(0),2)*
            (-1 + sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                       pow(v1(2) - v2(2),2))/l0))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) +
           (2*ks*pow(v1(0) - v2(0),2))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2)) +
           (2*ks*l0*(-1 + sqrt(pow(v1(0) - v2(0),2) +
                               pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2))/l0))/
           sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                pow(v1(2) - v2(2),2));
}
double ddEs_dx2dy1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (2*ks*l0*(v1(0) - v2(0))*(v1(1) - v2(1))*
            (-1 + sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                       pow(v1(2) - v2(2),2))/l0))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) -
           (2*ks*(v1(0) - v2(0))*(v1(1) - v2(1)))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2));
}
double ddEs_dx2dy2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (-2*ks*l0*(v1(0) - v2(0))*(v1(1) - v2(1))*
            (-1 + sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                       pow(v1(2) - v2(2),2))/l0))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) +
           (2*ks*(v1(0) - v2(0))*(v1(1) - v2(1)))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2));
}
double ddEs_dx2dz1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (2*ks*l0*(v1(0) - v2(0))*(-1 +
                                     sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                                          pow(v1(2) - v2(2),2))/l0)*(v1(2) - v2(2)))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) -
           (2*ks*(v1(0) - v2(0))*(v1(2) - v2(2)))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2));
}
double ddEs_dx2dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (-2*ks*l0*(v1(0) - v2(0))*(-1 +
                                      sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                                           pow(v1(2) - v2(2),2))/l0)*(v1(2) - v2(2)))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) +
           (2*ks*(v1(0) - v2(0))*(v1(2) - v2(2)))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2));
}

double ddEs_dy1dy1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (-2*ks*l0*pow(v1(1) - v2(1),2)*
            (-1 + sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                       pow(v1(2) - v2(2),2))/l0))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) +
           (2*ks*pow(v1(1) - v2(1),2))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2)) +
           (2*ks*l0*(-1 + sqrt(pow(v1(0) - v2(0),2) +
                               pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2))/l0))/
           sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                pow(v1(2) - v2(2),2));
}
double ddEs_dy1dy2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (2*ks*l0*pow(v1(1) - v2(1),2)*
            (-1 + sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                       pow(v1(2) - v2(2),2))/l0))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) -
           (2*ks*pow(v1(1) - v2(1),2))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2)) -
           (2*ks*l0*(-1 + sqrt(pow(v1(0) - v2(0),2) +
                               pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2))/l0))/
           sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                pow(v1(2) - v2(2),2));
}
double ddEs_dy1dz1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (-2*ks*l0*(v1(1) - v2(1))*(-1 +
                                      sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                                           pow(v1(2) - v2(2),2))/l0)*(v1(2) - v2(2)))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) +
           (2*ks*(v1(1) - v2(1))*(v1(2) - v2(2)))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2));
}
double ddEs_dy1dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (2*ks*l0*(v1(1) - v2(1))*(-1 +
                                     sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                                          pow(v1(2) - v2(2),2))/l0)*(v1(2) - v2(2)))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) -
           (2*ks*(v1(1) - v2(1))*(v1(2) - v2(2)))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2));
}

double ddEs_dy2dy2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (-2*ks*l0*pow(v1(1) - v2(1),2)*
            (-1 + sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                       pow(v1(2) - v2(2),2))/l0))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) +
           (2*ks*pow(v1(1) - v2(1),2))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2)) +
           (2*ks*l0*(-1 + sqrt(pow(v1(0) - v2(0),2) +
                               pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2))/l0))/
           sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                pow(v1(2) - v2(2),2));
}
double ddEs_dy2dz1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (2*ks*l0*(v1(1) - v2(1))*(-1 +
                                     sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                                          pow(v1(2) - v2(2),2))/l0)*(v1(2) - v2(2)))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) -
           (2*ks*(v1(1) - v2(1))*(v1(2) - v2(2)))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2));
}
double ddEs_dy2dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (-2*ks*l0*(v1(1) - v2(1))*(-1 +
                                      sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                                           pow(v1(2) - v2(2),2))/l0)*(v1(2) - v2(2)))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) +
           (2*ks*(v1(1) - v2(1))*(v1(2) - v2(2)))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2));
}

double ddEs_dz1dz1(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (2*ks*l0*(-1 + sqrt(pow(v1(0) - v2(0),2) +
                               pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2))/l0))/
           sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                pow(v1(2) - v2(2),2)) -
           (2*ks*l0*(-1 + sqrt(pow(v1(0) - v2(0),2) +
                               pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2))/l0)*
            pow(v1(2) - v2(2),2))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) +
           (2*ks*pow(v1(2) - v2(2),2))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2));
}
double ddEs_dz1dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (-2*ks*l0*(-1 + sqrt(pow(v1(0) - v2(0),2) +
                                pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2))/l0))/
           sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                pow(v1(2) - v2(2),2)) +
           (2*ks*l0*(-1 + sqrt(pow(v1(0) - v2(0),2) +
                               pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2))/l0)*
            pow(v1(2) - v2(2),2))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) -
           (2*ks*pow(v1(2) - v2(2),2))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2));
}
double ddEs_dz2dz2(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const double& l0, const double& ks) {
    return (2*ks*l0*(-1 + sqrt(pow(v1(0) - v2(0),2) +
                               pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2))/l0))/
           sqrt(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
                pow(v1(2) - v2(2),2)) -
           (2*ks*l0*(-1 + sqrt(pow(v1(0) - v2(0),2) +
                               pow(v1(1) - v2(1),2) + pow(v1(2) - v2(2),2))/l0)*
            pow(v1(2) - v2(2),2))/
           pow(pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
               pow(v1(2) - v2(2),2),1.5) +
           (2*ks*pow(v1(2) - v2(2),2))/
           (pow(v1(0) - v2(0),2) + pow(v1(1) - v2(1),2) +
            pow(v1(2) - v2(2),2));
}

/*
 *  ==========================================
 *      Shear Energy, 1st derivatives
 *  ==========================================
 */



/*
 *  ==========================================
 *      Shear Energy, 2nd derivatives
 *  ==========================================
 */



/*
 *  ==========================================
 *      Bending Energy, 1st derivatives
 *  ==========================================
 */



/*
 *  ==========================================
 *      Bending Energy, 2nd derivatives
 *  ==========================================
 */
>>>>>>> 760a66da6dd407a9409908295fcce13f5c25e423
