#ifndef PLATES_SHELLS_PARAMETERS_H
#define PLATES_SHELLS_PARAMETERS_H

class Parameters {
public:
    Parameters();

    // accessor
    int           ndof() const;
    int           nen() const;
    int           nst() const;
    int           iter_lim() const;
    double        dt() const;
    double        rec_len() const;
    double        rec_wid() const;
    unsigned long num_nodes_len() const;
    unsigned long num_nodes_wid() const;
    unsigned long nn() const;
    unsigned long nel() const;
    double        E_modulus() const;
    double        nu() const;
    double        rho() const;
    double        thk() const;
    double        vis() const;
    double        kstretch() const;
    double        kshear() const;
    double        kbend() const;

    // modifier
    void set_ndof(const int& var);
    void set_nen(const int& var);
    void set_nst(const int& var);
    void set_iter_lim(const int& var);
    void set_dt(const double& var);
    void set_rec_len(const double& var);
    void set_rec_wid(const double& var);
    void set_num_nodes_len(const unsigned long& var);
    void set_num_nodes_wid(const unsigned long& var);
    void set_nn(const unsigned long& var);
    void set_nel(const unsigned long& var);
    void set_E_modulus(const double& var);
    void set_nu(const double& var);
    void set_rho(const double& var);
    void set_thk(const double& var);
    void set_vis(const double& var);
    void set_kstretch();
    void set_kshear();
    void set_kbend();

    // other functions
    void print_parameters();

private:
    int             ndof_;                       // degree of freedom
    int             nen_;                        // number of nodes per element
    int             nst_;                        // total number of steps
    int             iter_lim_;                   // maximum number of iterations allowed per time step
    double          dt_;                         // step size
    double          rec_len_;                    // length of the rectangular domain
    double          rec_wid_;                    // width of the rectangular domain
    unsigned long   num_nodes_len_;              // number of nodes along the length
    unsigned long   num_nodes_wid_;              // number of nodes along the width
    unsigned long   nn_;                         // total number of nodes
    unsigned long   nel_;                        // total number of elements
    double          E_modulus_;                  // Young's modulus
    double          nu_;                         // Poisson's ratio
    double          rho_;                        // density
    double          thk_;                        // thickness
    double          vis_;                        // viscosity
    double          ks_;                         // stretch stiffness
    double          ksh_;                        // shearing stiffness
    double          kb_;                         // bending stiffness
};

#endif //PLATES_SHELLS_PARAMETERS_H
