//
//  fn_orthopoly.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 12.01.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_fn_orthopoly_hpp
#define UZLMathLib_fn_orthopoly_hpp

UZLMATH_NAMESPACE(orthoPoly)

// Legendre polynomial for x of degree n
double legendre(const int& n, const double& x);

// Associated Legendre polynomial
double assoc_legendre(const int& l, const int& m, const double& x);

UZLMATH_NAMESPACE_END

#endif
