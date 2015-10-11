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

/*!
 * @brief       The Spherical harmonics \f$Y^m_l(\theta, \phi)\f$.
 * @details     The Spherical harmonics \f$Y^m_l(\theta)\f$ for 
 *              \f$\theta\in[0,\pi]\f$, \f$\phi\in[0,2\pi)\f$, 
 *              \f$|m|\leq l\f$ is defined as
 *              \f{eqnarray*}{
 *                  Y^m_l(\theta, \phi) := \sqrt{\frac{2l + 1}{4\pi}\cdot\frac{(l-m)!}{(l+m)!}}\cdot P^m_l(\cos\theta)\cdot e^{\mathrm i m\phi},
 *              \f}
 *              where \f$P^m_l(x)\f$ denotes the associated legendre
 *              polynomial.
 *
 * @param[in]   m Degree of the spherical harmonic
 * @param[in]   l Order of the spherical harmonic
 * @param[in]   theta Colatitudinal coordinate
 * @param[in]   phi Longitudinal coordinate
 *
 * @return      The value of the spherical harmonic
 *
 * @sa          assoc_legendre
 *
 * @author      Denis-Michael Lux <denis.lux@icloud.com>
 * @date        11.10.15
 *
 * @since       0.1.1
 */
complex< double > Ylm(const int& l, const int& m, const double& theta, const double& phi)
{
    // norm factor
    double norm = (2.0 * l + 1.0) / (4.0 * constants< double >::pi);
    
    // needed indices
    int i;
    
    // evaluate factorials. First divide the sign in
    // each iteration and the multiply each factorial
    // step to the norm factor to prevent bitoverflow.
    // At last compute square root.
    for (i = 1; i <= l + m; ++i)
    {
        norm /= i;
    }
    
    for (i = 1; i <= l - m; ++i)
    {
        norm *= i;
    }
    
    norm = sqrt(norm);
    
    // complex term
    complex< double > e;
    e.polar(1, m * phi);
    
    // return value by multiplying the norm with the legendre
    // polynomial and the complex term
    return complex< double >(norm * orthoPoly::assoc_legendre(l, m, cos(theta)), 0) * e;
}

UZLMATH_NAMESPACE_END

#endif