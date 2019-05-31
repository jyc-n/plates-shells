#ifndef PLATES_SHELLS_LOADBC_H
#define PLATES_SHELLS_LOADBC_H

#include <vector>
#include <Eigen/Dense>

class Parameters;
class Geometry;

// TODO: add nonhomogeneous Dirichlet dofs options
struct Sets {
    Sets(bool type)
        : m_type(type) ,m_fixed({false, false, false})
    {}
    void setFixed(bool x, bool y, bool z) {
        m_fixed[0] = x;
        m_fixed[1] = y;
        m_fixed[2] = z;
    }
    bool m_type;                // true - force sets, false - displacement sets
    std::vector<bool> m_fixed;  // true - fixed in 0-x, 1-y, 2-z
    std::vector<int>  m_nodes;
};

class Boundary {

public:
    Boundary(Parameters* SimPar, Geometry* SimGeo);
    ~Boundary();

    // initialize boundary conditions
    void initBC();

    // check if a given dof with given index is a Dirichlet dof
    bool inDirichletBC(int index);

    // public variables that will be used by solver
    std::vector<int>  m_dirichletDofs;               // list of dofs in Dirichlet BC
    Eigen::VectorXd m_fext;                          // external force vector

private:
    void configForce(const Sets& t_set);
    void configGravity();
    void findDirichletDofs(const Sets& t_set);

    Parameters* m_SimPar;
    Geometry*   m_SimGeo;

    bool ENABLE_GRAVITY;
};

#endif //PLATES_SHELLS_LOADBC_H
