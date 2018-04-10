#include <iostream>
#include "parameters.h"

// default constructor
Parameters::Parameters() {
    ndof_               =   0;
    nen_                =   0;
    rec_len_            =   0;
    rec_wid_            =   0;
    num_nodes_len_      =   0;
    num_nodes_wid_      =   0;
    nn_                 =   0;
    nel_                =   0;
}

// accessors
int           Parameters::ndof() const          { return ndof_; }
int           Parameters::nen() const           { return nen_; }
double        Parameters::rec_len() const       { return rec_len_; }
double        Parameters::rec_wid() const       { return rec_wid_; }
unsigned long Parameters::num_nodes_len() const { return num_nodes_len_; }
unsigned long Parameters::num_nodes_wid() const { return num_nodes_wid_; }
unsigned long Parameters::nn() const            { return nn_; }
unsigned long Parameters::nel() const           { return nel_; }

// modifiers
void Parameters::set_ndof(const int &var)                    { ndof_ = var; }
void Parameters::set_nen(const int &var)                     { nen_ = var; }
void Parameters::set_rec_len(const double &var)              { rec_len_ = var; }
void Parameters::set_rec_wid(const double &var)              { rec_wid_ = var; }
void Parameters::set_num_nodes_len(const unsigned long &var) { num_nodes_len_ = var; }
void Parameters::set_num_nodes_wid(const unsigned long &var) { num_nodes_wid_ = var; }
void Parameters::set_nn(const unsigned long &var)            { nn_ = var; }
void Parameters::set_nel(const unsigned long &var)           { nel_ = var; }

void Parameters::print_parameters() {
    std::cout << "List of parameters" << '\n';
    std::cout << "length of the rectangular domain: " << rec_len_ << '\n';
    std::cout << "width of the rectangular domain: " << rec_wid_ << '\n';
    std::cout << "number of nodes along the length: " << num_nodes_len_ << '\n';
    std::cout << "number of nodes along the width: " << num_nodes_wid_ << '\n';
}

