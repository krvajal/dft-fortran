//
// Created by Miguel Angel Carvajal on 8/31/16.
//

#ifndef DFT_QUAD_SIMPSON_H_H
#define DFT_QUAD_SIMPSON_H_H
#include <armadillo>

/// Perform a cuadrature calculation of the function using simpsom method
/// \param v  samples of the function to be integrated
/// \return the value of the integral of the function

template <typename T>
T simpson(const arma::Col<T> &v){

    int size = v.n_rows;
    double sum = v(0);
    for(int i = 1; i < size -1; i++){
        int fact = 1 << (i%2 + 1);
        sum += fact * v(i);

    }
    sum += v(size -1);
    return sum/3.;

}


#endif //DFT_QUAD_SIMPSON_H_H
