#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <Eigen/Dense>
#include "solver.h"
#include "element.h"
#include "geometry.h"

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

void time_solver(Parameters& SimPar) {
    //m_node_lst = new Node[SimPar.nn()];
    Node m_node_lst[SimPar.nn()];
    Element m_element_lst[SimPar.nel()];

    prepare_solver(SimPar, m_node_lst, m_element_lst);
    // solver codes goes here

}

// initialize a list for
void prepare_solver(Parameters& SimPar, Node* l_node, Element* l_element) {
    for (int i = 0; i < SimPar.nn(); i++) {
        Node temp(m_coord(i,0), m_coord(i,1), m_coord(i,2));
        l_node[i] = temp;
    }

    for (int i = 0; i < SimPar.nel(); i++) {
        Element temp(i+1, m_conn(i,0), m_conn(i,1), m_conn(i,2));

        // pointers to 3 node object
        Node* pn1 = &l_node[m_conn(i,0)-1];
        Node* pn2 = &l_node[m_conn(i,1)-1];
        Node* pn3 = &l_node[m_conn(i,2)-1];
        temp.set_node(pn1, pn2, pn3);

        temp.calculate_dir();
        temp.calculate_normal();
        temp.calculate_angle();
        temp.find_nearby_element(SimPar);
        l_element[i] = temp;
    }
}