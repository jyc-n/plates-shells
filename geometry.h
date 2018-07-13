#ifndef PLATES_SHELLS_GEOMETRY_H
#define PLATES_SHELLS_GEOMETRY_H

#include <Eigen/Dense>

class Parameters;
class Node;
class Element;

class Geometry {

public:
    Geometry(Parameters* SimPar);
    ~Geometry();

    // modifier
    void set_ndof(const int var);
    void set_nen(const int var);
    void set_rec_len(const double var);
    void set_rec_wid(const double var);
    void set_num_nodes_len(const unsigned int var);
    void set_num_nodes_wid(const unsigned int var);
    void set_nn(const unsigned int var);
    void set_nel(const unsigned int var);

    void buildGeo();
    void calcMass();

    // accessor
    int           ndof() const;
    int           nen() const;
    double        rec_len() const;
    double        rec_wid() const;
    unsigned int  num_nodes_len() const;
    unsigned int  num_nodes_wid() const;
    unsigned int  nn() const;
    unsigned int  nel() const;

    void printGeo();

    // public geometry variables
    Eigen::MatrixXi m_conn;               // connectivity matrix     (nel x nen)
    Eigen::MatrixXd m_coord;              // coordinates matrix      (nn  x nsd)
    Eigen::VectorXd m_dof;                // degree of freedom x1, y1, (z1), x2, y2, (z2), ... as a vector
    Eigen::MatrixXi map_nodes;            // map of nodes, for verification purpose
    Eigen::VectorXd m_mass;               // nodal mass vector, mx1, my1, mz1, mx2, my2, mz2, ...

    Node*    m_nodeList;
    Element* m_elementList;

private:

    // geometry parameters
    int             m_ndof;                       // degree of freedom
    int             m_nen;                        // number of nodes per element
    double          m_rec_len;                    // length of the rectangular domain
    double          m_rec_wid;                    // width of the rectangular domain
    unsigned int    m_num_nodes_len;              // number of nodes along the length
    unsigned int    m_num_nodes_wid;              // number of nodes along the width
    unsigned int    m_nn;                         // total number of nodes
    unsigned int    m_nel;                        // total number of elements

    // pointer to parameters
    Parameters* m_SimPar;

};

// inline implementations for short functions
inline void Geometry::set_ndof(const int var)                    { m_ndof = var; }
inline void Geometry::set_nen(const int var)                     { m_nen = var; }
inline void Geometry::set_rec_len(const double var)              { m_rec_len = var; }
inline void Geometry::set_rec_wid(const double var)              { m_rec_wid = var; }
inline void Geometry::set_num_nodes_len(const unsigned int var)  { m_num_nodes_len = var; }
inline void Geometry::set_num_nodes_wid(const unsigned int var)  { m_num_nodes_wid = var; }
inline void Geometry::set_nn(const unsigned int var)             { m_nn = var; }
inline void Geometry::set_nel(const unsigned int var)            { m_nel = var; }

inline int           Geometry::ndof() const                      { return m_ndof; }
inline int           Geometry::nen() const                       { return m_nen; }
inline double        Geometry::rec_len() const                   { return m_rec_len; }
inline double        Geometry::rec_wid() const                   { return m_rec_wid; }
inline unsigned int  Geometry::num_nodes_len() const             { return m_num_nodes_len; }
inline unsigned int  Geometry::num_nodes_wid() const             { return m_num_nodes_wid; }
inline unsigned int  Geometry::nn() const                        { return m_nn; }
inline unsigned int  Geometry::nel() const                       { return m_nel; }

#endif //PLATES_SHELLS_GEOMETRY_H
