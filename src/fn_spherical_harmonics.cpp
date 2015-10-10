//
//  fn_spherical_harmonics.cpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 10.10.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_fn_spherical_harmonics_cpp
#define UZLMathLib_fn_spherical_harmonics_cpp

#include <uzlmath>

UZLMATH_NAMESPACE(spharmonics)

complex< double > Ylm(const int& l, const int& m, const double& nu, const double& phi)
{
    // scale factor
    double scale = 1 / (2 * constants< double >::pi);
    
    // compute factorials
    double lmm = 1.0;
    double lpm = 1.0;
    
    // used indices
    int i;
    
    // compute factorials
    for (i = 2; i <= l - m; ++i)
    {
        lpm *= i;
        lmm *= i;
    }
    
    for (i = 1; i <= 2 * m; ++i)
    {
        lpm *= (l - m + i);
    }
    
    // compute norm factor
    double norm = sqrt( (2.0 * l + 1.0) * lmm / lpm / 2.0 );
    
    // get complex value
    complex< double > e;
    e.polar(1, static_cast< double >(m) * phi);
    
    std::cout << " IN SPHARMONICS " << e << std::endl;
    
    return complex< double >(scale * norm * orthoPoly::assoc_legendre(l, m, cos(nu)), 0) * e;
}

UZLMATH_NAMESPACE_END

#endif