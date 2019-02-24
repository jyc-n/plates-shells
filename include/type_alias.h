#ifndef TYPE_ALIAS_H
#define TYPE_ALIAS_H

#define EIGEN_USE_NEW_STDVECTOR
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/StdVector>

using VectorNodes = std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> >;
using VectorMesh = std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i> >;
using VectorN = Eigen::VectorXd;
using SpMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using SparseEntries = std::vector<Eigen::Triplet<double> >;

#endif // TYPE_ALIAS_H