#ifndef PLATES_SHELLS_LOADBC_H
#define PLATES_SHELLS_LOADBC_H

#include <vector>
#include <Eigen/Dense>

class Parameters;
class Geometry;

struct Sets {
    bool m_free = true;         // false - BC applied to this set, true - no BC applied
    int m_typeBC = 0;           // 0 - free, 1 - pinned in x, 2 - pinned in y, 3 - clamped, 4 - displacement applied (not implemented yet)
    double m_force[3] = {0};    // force in x,y,z direction
    std::vector<int> m_nodes;
};

class Boundary {

public:
    Boundary(Parameters* SimPar, Geometry* SimGeo);
    ~Boundary();

    void initBC();
    void getSets();
    void configBC();
    void fixDir(int num);
    void buildBCinfo();
    void getSpecifiedDof();

    // public variables that will be used by solver
    int m_numTotal;
    int m_numFree;
    std::vector<int>  m_specifiedDof;               // list of specified dof
    Eigen::VectorXd m_fext;                         // external force vector

private:
    Parameters* m_SimPar;
    Geometry*   m_SimGeo;

    int m_bounds;
    bool m_2D;

    std::vector<Sets*> m_sets;                      // node sets of different boundaries
    std::vector<bool> m_disp;                       // true - free dofs, false - specified dofs
};

#endif //PLATES_SHELLS_LOADBC_H
