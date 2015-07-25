//
//  fn_dwt.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 03.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_fn_dwt_hpp
#define UZLMathLib_fn_dwt_hpp

UZLMATH_NAMESPACE(DWT)
    
/*- For more information/implementation details see fn_dwt.cpp file! -*/

// Quadrature weights for the discrete Wigner transform
auto quadrature_weights(const int& bandwidth) -> vector< double >;

// Weighted Wigner-d (L^2 normalized Wigner-d entries) matrix for the discrete Wigner transform
auto weighted_wigner_d_matrix(const int& bandwidth, const int& M, const int& Mp, const vector< double >& weights) -> matrix< double >;

// Wigner-d (L^2 normalized Wigner-d entries) matrix for the discrete Wigner transform
auto wigner_d_matrix(const int& bandwidth, const int& M, const int& Mp) -> matrix< double >;

UZLMATH_NAMESPACE_END

#endif
