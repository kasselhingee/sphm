#include "SvMF_links.h"

mata1 alignedGcpp(veca1 m, mata1 P) {
    int p = m.size();
    
    // Compute m %*% t(m) which is the outer product of m with itself
    mata1 mproj = m * m.transpose();

    // first remove the mean direction from all directions in P. P 'no m'
    mata1 Pnom = (Eigen::MatrixXd::Identity(p, p) - mproj) * P;

    // Remove the mean direction from all directions in P
    for (int j = 1; j < p - 1; ++j) {
        // Normalize Pnom[, j]
        Pnom.col(j) /= Pnom.col(j).norm();

        // Update Pnom by removing the jth direction from subsequent columns
        Pnom.block(0, j + 1, p, p - j - 1) -= Pnom.col(j) * Pnom.col(j).transpose() * Pnom.block(0, j + 1, p, p - j - 1);
    }

    // Normalize the last column of Pnom since it misses out in the above loop
    Pnom.col(p - 1) /= Pnom.col(p - 1).norm();

    // Replace the first column of Pnom with m
    Pnom.col(0) = m;

    return Pnom;
}
