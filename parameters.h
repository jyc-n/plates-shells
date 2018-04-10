#ifndef PLATES_SHELLS_PARAMETERS_H
#define PLATES_SHELLS_PARAMETERS_H

class Parameters {
public:
    Parameters();

    // accessor
    int           ndof() const;
    int           nen() const;
    double        rec_len() const;
    double        rec_wid() const;
    unsigned long num_nodes_len() const;
    unsigned long num_nodes_wid() const;
    unsigned long nn() const;
    unsigned long nel() const;

    // modifier
    void set_ndof(const int& var);
    void set_nen(const int& var) ;
    void set_rec_len(const double& var);
    void set_rec_wid(const double& var);
    void set_num_nodes_len(const unsigned long& var);
    void set_num_nodes_wid(const unsigned long& var);
    void set_nn(const unsigned long& var);
    void set_nel(const unsigned long& var);

    // other functions
    void print_parameters();

private:
    int             ndof_;                       // degree of freedom
    int             nen_;                        // number of nodes per element
    double          rec_len_;                    // length of the rectangular domain
    double          rec_wid_;                    // width of the rectangular domain
    unsigned long   num_nodes_len_;              // number of nodes along the length
    unsigned long   num_nodes_wid_;              // number of nodes along the width
    unsigned long   nn_;                         // total number of nodes
    unsigned long   nel_;                        // total number of elements
};

#endif //PLATES_SHELLS_PARAMETERS_H
