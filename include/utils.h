//
// Created by jakob on 4/13/25.
//

#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>


namespace bisam {
    arma::ivec repeat_value(int value, int times);

    arma::mat kronecker_product(const arma::mat &A, const arma::mat &B);
}
#endif //UTILS_H
