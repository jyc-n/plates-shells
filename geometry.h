#ifndef PLATES_SHELLS_GEOMETRY_H
#define PLATES_SHELLS_GEOMETRY_H

#include <Eigen/Dense>

class Node;
class Element;

// global geometry variables
extern Eigen::MatrixXi m_conn;               // connectivity matrix     (nel x nen)
extern Eigen::MatrixXd m_coord;              // coordinates matrix      (nn  x nsd)
extern Eigen::VectorXd m_dof;                // degree of freedom x1, y1, (z1), x2, y2, (z2), ... as a vector
extern Eigen::MatrixXi map_nodes;            // map of nodes, for verification purpose

#endif //PLATES_SHELLS_GEOMETRY_H
