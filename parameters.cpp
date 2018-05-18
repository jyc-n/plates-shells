#include <iostream>
#include <cmath>
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
int           Parameters::nst() const           { return nst_; }
int           Parameters::iter_lim() const      { return iter_lim_; }
double        Parameters::dt() const            { return dt_; }
double        Parameters::rec_len() const       { return rec_len_; }
double        Parameters::rec_wid() const       { return rec_wid_; }
unsigned long Parameters::num_nodes_len() const { return num_nodes_len_; }
unsigned long Parameters::num_nodes_wid() const { return num_nodes_wid_; }
unsigned long Parameters::nn() const            { return nn_; }
unsigned long Parameters::nel() const           { return nel_; }
double        Parameters::E_modulus() const     { return E_modulus_; }
double        Parameters::nu() const            { return nu_; }
double        Parameters::rho() const           { return rho_; }
double        Parameters::thk() const           { return thk_; }
double        Parameters::vis() const           { return vis_; }
double        Parameters::kstretch() const      { return ks_; }
double        Parameters::kshear() const        { return ksh_; }
double        Parameters::kbend() const         { return kb_; }

// modifiers
void Parameters::set_ndof(const int &var)                    { ndof_ = var; }
void Parameters::set_nen(const int &var)                     { nen_ = var; }
void Parameters::set_nst(const int &var)                     { nst_ = var; }
void Parameters::set_iter_lim(const int &var)                { iter_lim_ = var; }
void Parameters::set_dt(const double &var)                   { dt_ = var; }
void Parameters::set_rec_len(const double &var)              { rec_len_ = var; }
void Parameters::set_rec_wid(const double &var)              { rec_wid_ = var; }
void Parameters::set_num_nodes_len(const unsigned long &var) { num_nodes_len_ = var; }
void Parameters::set_num_nodes_wid(const unsigned long &var) { num_nodes_wid_ = var; }
void Parameters::set_nn(const unsigned long &var)            { nn_ = var; }
void Parameters::set_nel(const unsigned long &var)           { nel_ = var; }
void Parameters::set_E_modulus(const double& var)            { E_modulus_ = var; }
void Parameters::set_nu(const double& var)                   { nu_ = var; }
void Parameters::set_rho(const double& var)                  { rho_ = var; }
void Parameters::set_thk(const double& var)                  { thk_ = var; }
void Parameters::set_vis(const double& var)                  { vis_ = var; }
void Parameters::set_kstretch()                              { ks_ = E_modulus_ * thk_; }
void Parameters::set_kshear()                                { ksh_ = E_modulus_ * thk_; }
void Parameters::set_kbend()                                 { kb_ = E_modulus_ * pow(thk_,3) / (12.0 * (1 - pow(nu_,2))); }

void Parameters::print_parameters() {
    std::cout << "List of parameters" << '\n';
    std::cout << "# of Nodes:       " << nn_ << '\n';
    std::cout << "# of Elements:    " << nel_ << '\n';
    std::cout << "Step size:        " << dt_ << '\n';
    std::cout << "Total # of steps: " << nst_ << '\n';
    std::cout << "---------------------------------------" << std::endl;
}