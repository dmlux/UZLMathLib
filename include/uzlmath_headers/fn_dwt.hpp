//
//  fn_dwt.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 03.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_fn_dwt_hpp
#define uzlmath_fn_dwt_hpp

/*! 
 * @namespace   dwt
 * @brief       A namespace that contains all functions to generate matrices and
 *              weights for the Discrete Wigner Transform (_DWT_) of a data vector
 *              \f$s\f$
 * @details     The functions that are described in this namespace are mainly used
 *              to construct the matrices that are used for the DWT that is described
 *              in the paper 'FFTs on the Rotation Group' by P. J. Kostelec and D. N.
 *              Rockmore. The DWT is defined as a collection of sums of the form
 *              \f[
 *                  \hat{s}(l,M,M')=\sum\limits_{k=0}^{2B-1}w_B(k)\tilde{d}^l_{M,M'}(\beta_k)[s]_k
 *              \f]
 *              where \f$\max(|M|,|M'|)\leq l\leq B\f$, \f$\tilde{d}^l_{M,M'}\f$ is
 *              a normalized Wigner d-function of degree \f$l\f$ and orders \f$M\f$,
 *              \f$M'\f$ and \f$\beta_k = \frac{\pi(2k+1)}{4B}\f$.
 *
 * @since       0.0.1
 *
 * @author      Denis-Michael Lux <denis.lux@icloud.com>
 * @date        03.05.15
 */
namespace DWT
{
    
    // For more information/implementation details see fn_dwt.cpp file!
    
    // Quadrature weights for the discrete Wigner transform
    auto quadrature_weights(const unsigned int& bandwidth) -> vector< double >;
    
    // Weighted Wigner-d (L^2 normalized Wigner-d entries) matrix for the discrete Wigner transform
    auto weighted_wigner_d_matrix(const unsigned int& bandwidth, const int& M, const int& Mp, const vector< double >& weights) -> matrix< double >;
    
    // Wigner-d (L^2 normalized Wigner-d entries) matrix for the discrete Wigner transform
    auto wigner_d_matrix(const unsigned int& bandwidth, const int& M, const int& Mp) -> matrix< double >;
    
}

#endif
