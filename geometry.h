#ifndef PLATES_SHELLS_GEOMETRY_H
#define PLATES_SHELLS_GEOMETRY_H

#include <Eigen/Dense>

#include "parameters.h"

class Geometry {
public:

    // modifier
    void init_lst_nodes(Parameters& sims);
    void init_lst_elements(Parameters& sims);

    // accessor
    Eigen::MatrixXd lst_coord() const;            // return coordinates table
    Eigen::MatrixXi lst_conn() const;             // return connectivity table

    void print_geo(Parameters& sims);

private:

    Eigen::MatrixXi lst_elements;           // connectivity table
    Eigen::MatrixXd lst_nodes;              // coordinates table
    Eigen::MatrixXi map_nodes;              // map of nodes, for verification purpose
};

#endif //PLATES_SHELLS_GEOMETRY_H
