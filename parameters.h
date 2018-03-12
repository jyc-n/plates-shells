#ifndef PLATES_SHELLS_PARAMETERS_H
#define PLATES_SHELLS_PARAMETERS_H

#include <iostream>
#include <fstream>
#include <string>

class Parameters {
public:
    Parameters();

    // read input parameters
    void print_parameters();

    // modifiers
    void set_domain_len(const double& len) { rec_len_ = len; };
    void set_domain_wid(const double& wid) { rec_wid_ = wid; };
    void set_num_nodes_len(const unsigned int& n_len) { num_nodes_len_ = n_len; };
    void set_num_nodes_wid(const unsigned int& n_wid) { num_nodes_wid_ = n_wid; };
    void set_dof(const int& dof) { dof_ = dof; }
    void set_nen(const int& nen) { nen_ = nen; }

    // accessors
    double& rec_len() { return rec_len_; }
    double& rec_wid() { return rec_wid_; }
    unsigned int& num_node_len() { return num_nodes_len_; }
    unsigned int& num_node_wid() { return num_nodes_wid_; }
    int& get_dof() { return dof_; }
    int& get_nen() { return nen_; }

private:
    int dof_, nen_;
    double rec_len_, rec_wid_;
    unsigned int num_nodes_len_, num_nodes_wid_;
};

#endif //PLATES_SHELLS_PARAMETERS_H
