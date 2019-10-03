# plates-shells

## Description

This is a numerical simulation solver for thin-flexible structures (such as cloth, plates, shells). The code was written primarily based on the following papers:

1. Grinspun, E., Hirani, A. N., Desbrun, M., & Schroeder, P. (2003). Discrete Shells. *Proceedings of the 2003 ACM SIGGRAPH*, 62–67. Retrieved from http://dl.acm.org/citation.cfm?id=846276.846284
2. Tamstorf, R., & Grinspun, E. (2013). Discrete bending forces and their Jacobians. *Graphical Models*, *75*(6), 362–370. https://doi.org/10.1016/j.gmod.2013.07.001



The main difference between each branch is the way that the elastic force and jacobian were derived.

#### Master Branch 
- Derivatives were analytically derived
- Modifications should be made onto this branch

**The following branches are deprecated**

i.e. don't use them

**Symbolic Branch (deprecated)**
Energy derivatives calculated using symbolic math tools

**Numerical Branch (deprecated)**
Numerical differentiation (central difference) were used

## Set Up

To use this solver, several **external libraries** are needed:

1. [CMake](https://cmake.org/download/)

   Any version of CMake can be used. Latest version is preferred.

   Also, if the `g++/clang++` compiler is new enough to be compatible with `C++17`, automatic I/O path benefited from the `filesystem` library can be used. Modify the `arguments.h` to do this.

2. [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) 

   After downloading **Eigen**, unzip the file, rename the fold (for instance, Eigen), and move the entire folder anywhere you like. Create the soft link in your system include direction. 

   For example, your Eigen (version 3.3.7) folder is in

   `/usr/local/Cellar`

   Create a soft link (works for Mac/Linux)

   `ln -s /usr/local/Cellar/Eigen/3.3.7/include/eigen3/Eigen /usr/local/include/Eigen`

3. [Pardiso](https://www.pardiso-project.org) 

    (refer to their website for set up instructions)

