#include <iostream>
#include <cmath>
#include "parameters.h"
#include "geometry.h"

// default constructor
Parameters::Parameters() {
    m_SimGeo = nullptr;
    m_inputPath = "/Users/chenjingyu/Dropbox/Research/Codes/plates-shells/";
    m_outputPath = "/Users/chenjingyu/Dropbox/Research/Codes/plates-shells/results/";
}


// modifiers
void Parameters::link_geo(Geometry* geo)                    { m_SimGeo = geo; }
void Parameters::set_outop(const bool var)                  { outop_ = var; }
void Parameters::set_solver_op(const bool var)              { solver_op_ = var; }
void Parameters::set_info_style(const bool var)             { info_style_ = var; }
void Parameters::set_out_freq(const int var)                { out_freq_ = var; }
void Parameters::set_nst(const int var)                     { nst_ = var; }
void Parameters::set_iter_lim(const int var)                { iter_lim_ = var; }
void Parameters::set_dt(const double var)                   { dt_ = var; }
void Parameters::set_E_modulus(const double var)            { E_modulus_ = var; }
void Parameters::set_ctol(const double var)                 { ctol_ = var; }
void Parameters::set_nu(const double var)                   { nu_ = var; }
void Parameters::set_rho(const double var)                  { rho_ = var; }
void Parameters::set_thk(const double var)                  { thk_ = var; }
void Parameters::set_vis(const double var)                  { vis_ = var; }
void Parameters::set_gconst(const double var)               { gconst_ = var; }
void Parameters::set_kstretch()                             { ks_ = E_modulus_ * thk_; }
void Parameters::set_kshear()                               { ksh_ = E_modulus_ * thk_; }
void Parameters::set_kbend()                                { kb_ = E_modulus_ * pow(thk_,3) / (24.0 * (1 - pow(nu_,2))); }

// accessors
std::string   Parameters::inputPath()  const    { return m_inputPath; }
std::string   Parameters::outputPath() const    { return m_outputPath; }
bool          Parameters::outop() const         { return outop_; }
bool          Parameters::solver_op() const     { return solver_op_; }
int           Parameters::out_freq() const      { return out_freq_; }
int           Parameters::nst() const           { return nst_; }
int           Parameters::iter_lim() const      { return iter_lim_; }
double        Parameters::dt() const            { return dt_; }
double        Parameters::E_modulus() const     { return E_modulus_; }
double        Parameters::ctol() const          { return ctol_; }
double        Parameters::nu() const            { return nu_; }
double        Parameters::rho() const           { return rho_; }
double        Parameters::thk() const           { return thk_; }
double        Parameters::vis() const           { return vis_; }
double        Parameters::gconst() const        { return gconst_; }
double        Parameters::kstretch() const      { return ks_; }
double        Parameters::kshear() const        { return ksh_; }
double        Parameters::kbend() const         { return kb_; }

void Parameters::print_parameters() {
    std::cout << "\n---------------------------------------" << std::endl;
    std::cout << "\t\tList of parameters" << '\n';
    std::cout << "\t\t# of Nodes:       " << m_SimGeo->nn() << '\n';
    std::cout << "\t\t# of Elements:    " << m_SimGeo->nel() << '\n';
    std::cout << "\t\tStep size:        " << dt_ << '\n';
    std::cout << "\t\tTotal # of steps: " << nst_ << '\n';
    std::cout << "---------------------------------------\n" << std::endl;
}