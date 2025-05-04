#include "eigen3/Eigen/SparseCore"
#include "eigen3/Eigen/SparseLU"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/OrderingMethods"
#include <iostream>
#include "flint/flint.h"
#include "flint/fq_nmod.h"

typedef Eigen::SparseMatrix<int, Eigen::ColMajor> m_i;

int main(){
    // 1, 0, 1, 0, 0
    // 1, 1, 0, 0, 0
    // 0, 1, 0, 1, 0
    // 0, 0, 0, 6, 1
    m_i M{4, 5};

    M.insert(0, 0) = 1;
    M.insert(0, 2) = 1;
    M.insert(1, 0) = 1;
    M.insert(1, 1) = 1;
    M.insert(2, 1) = 1;
    M.insert(2, 3) = 1;
    M.insert(3, 3) = 6;
    M.insert(3, 4) = 1;

    std::cout << M << "\n";

    // auto A = M.block(0, 0, 3, 3);
    // std::cout << A;


    // Eigen::SparseLU<m_i, Eigen::COLAMDOrdering<int>> LU{M};

    // LU.analyzePattern(M);
    // LU.factorize(M);

    // auto U = LU.matrixU();
    // auto trueU = U.toSparse();


    
}