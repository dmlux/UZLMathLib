//
//  fn_spherical_harmonics.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 10.10.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_fn_spherical_harmonics_hpp
#define UZLMathLib_fn_spherical_harmonics_hpp

UZLMATH_NAMESPACE(spharmonics)

complex< double > Ylm(const int& l, const int& m, const double& nu, const double& phi);

UZLMATH_NAMESPACE_END

#endif
