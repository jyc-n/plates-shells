#ifndef PLATES_SHELLS_PARAMETERS_H
#define PLATES_SHELLS_PARAMETERS_H

#include <string>

class Geometry;

class Parameters {
public:
    Parameters();

    // modifier
    void link_geo(Geometry* geo);
    void set_outop(const bool var);
    void set_solver_op(const bool var);
    void set_info_style(const bool var);
    void set_out_freq(const int var);
    void set_nst(const int var);
    void set_iter_lim(const int var);
    void set_dt(const double var);
    void set_E_modulus(const double var);
    void set_ctol(const double var);
    void set_nu(const double var);
    void set_rho(const double var);
    void set_thk(const double var);
    void set_vis(const double var);
    void set_gconst(const double var);
    void set_kstretch();
    void set_kshear();
    void set_kbend();

    // accessor
    std::string   inputPath()  const;
    std::string   outputPath() const;
    bool          outop() const;
    bool          solver_op() const;
    bool          info_style() const;
    int           out_freq() const;
    int           nst() const;
    int           iter_lim() const;
    double        dt() const;
    double        E_modulus() const;
    double        ctol() const;
    double        nu() const;
    double        rho() const;
    double        thk() const;
    double        vis() const;
    double        gconst() const;
    double        kstretch() const;
    double        kshear() const;
    double        kbend() const;

    // other functions
    void print_parameters();

private:

    // pointer to the geometry class
    Geometry* m_SimGeo;

    // parameters
    bool            outop_;                      // output option
    bool            solver_op_;                  // solver option
    bool            info_style_;                 // console output style
    int             out_freq_;                   // output frequency
    int             nst_;                        // total number of steps
    int             iter_lim_;                   // maximum number of iterations allowed per time step
    double          dt_;                         // step size
    double          ctol_;                       // tolerance scaling function
    double          E_modulus_;                  // Young's modulus
    double          nu_;                         // Poisson's ratio
    double          rho_;                        // density
    double          thk_;                        // thickness
    double          vis_;                        // viscosity
    double          gconst_;
    double          ks_;                         // stretch stiffness
    double          ksh_;                        // shearing stiffness
    double          kb_;                         // bending stiffness

    std::string m_inputPath;
    std::string m_outputPath;
};

#endif //PLATES_SHELLS_PARAMETERS_H
