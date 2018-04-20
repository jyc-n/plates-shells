#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <Eigen/Dense>
#include "solver.h"
#include "element.h"

void fake_solver(Parameters& SimPar) {
    /*
    double ti = 0.0;
    double tf = 2.0;
    double dt = 0.1;
    double amp = 2.0 * dt / (tf - ti);

    Eigen::MatrixXd fake_motion = Eigen::MatrixXd::Zero(SimPar.nn(), SimPar.ndof());
    for (int i = 0; i < SimPar.nn(); i++) {
        for (int j = 0; j < SimPar.ndof(); j++) {
            if (j == 2) { // z-direction
                fake_motion(i,j) = amp * sin(M_PI / 4.0 * InitGeo.lst_coord()(i,0));
            }
        }
    }
    std::cout << fake_motion << '\n';

    Eigen::MatrixXd coord = InitGeo.lst_coord();
    int counter = 0;
    for (double t = ti; t < tf + 0.0000001; t += dt) {
        std::string filepath = "/Users/chenjingyu/git/plates-shells/";
        std::string filename = "result" + std::to_string(counter) + ".txt";
        std::ofstream myfile((filepath+filename).c_str());
        for (int i = 0; i < SimPar.nn(); i++) {
            for (int j = 0; j < SimPar.ndof(); j++) {
                myfile << std::setprecision(6) << std::fixed << coord(i,j) << '\t';
            }
            myfile << std::endl;
        }
        coord += fake_motion;
        counter++;
        //std::cout << '\n';
    }*/
}

void time_solver() {
    //Geometry InitGeo;
    //Parameters SimPar;
    //init_solver(InitGeo, SimPar);
}