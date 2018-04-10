#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <Eigen/Dense>
#include "solver.h"
#include "element.h"

void fake_solver(Geometry& InitGeo, Parameters& SimPar) {
    double ti = 0.0;
    double tf = 2.0;
    double dt = 0.4;
    double amp = 2.0 * dt / (tf - ti);

    Eigen::MatrixXd fake_motion = Eigen::MatrixXd::Zero(SimPar.nn(), SimPar.ndof());
    for (int i = 0; i < SimPar.nn(); i++) {
        for (int j = 0; j < SimPar.ndof(); j++) {
            fake_motion(i,j) = 2.0 * sin(M_PI / 2.0 * InitGeo.lst_coord()(i,0));
        }
    }
    //std::cout << fake_motion << '\n';

    Eigen::MatrixXd coord = InitGeo.lst_coord();
    int counter = 0;
    for (double t = ti; t < tf + 0.0000001; t += dt) {
        coord += fake_motion;

        counter++;
        std::string filepath = "/Users/chenjingyu/git/plates-shells/";
        std::string filename = "result" + std::to_string(counter) + ".txt";
        std::ofstream myfile((filepath+filename).c_str());
        for (int i = 0; i < SimPar.nn(); i++) {
            for (int j = 0; j < SimPar.ndof(); j++) {
                myfile << std::setprecision(6) << std::fixed << coord(i,j) << '\t';
            }
            myfile << std::endl;
        }
        //std::cout << '\n';
    }
}

void time_solver() {
    //Geometry InitGeo;
    //Parameters SimPar;
    //init_solver(InitGeo, SimPar);
}