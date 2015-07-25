//
//  fn_wigner_d.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 12.01.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_fn_wigner_d_hpp
#define UZLMathLib_fn_wigner_d_hpp

UZLMATH_NAMESPACE(wigner)

// Wigner-d function
auto wigner_d(const int& J, const int& M, const int& Mp, const double& beta) -> double;

// L^2 normalized Wiger-d function
auto wigner_d_l2normalized(const int& J, const int& M, const int& Mp, const double& beta) -> double;

// Wigner-D function
auto wigner_D(const int& J, const int& M, const int& Mp, const double& alpha, const double& beta, const double& gamma) -> const complex< double >;

// L^2 normalized Wigner-D function
auto wigner_D_l2normalized(const int& J, const int& M, const int& Mp, const double& alpha, const double& beta, const double& gamma) -> const complex< double >;

UZLMATH_NAMESPACE_END

#endif