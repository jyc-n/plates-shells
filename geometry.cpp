#include "geometry.h"
#include "parameters.h"

// global geometry variables
Eigen::MatrixXi m_conn;               // connectivity matrix     (nel x nen)
Eigen::MatrixXd m_coord;              // coordinates matrix      (nn  x nsd)
Eigen::VectorXd m_dof;                // degree of freedom x1, y1, (z1), x2, y2, (z2), ... as a vector
Eigen::MatrixXi map_nodes;            // map of nodes, for verification purpose