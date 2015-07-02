//
//  fn_cx_ntrt.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 12.01.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_fn_cx_ntrt_hpp
#define uzlmath_fn_cx_ntrt_hpp

/*!
 * @brief       Calculates the nth root of given number. And allows to calculate 
 *              root of negative real numbers.
 * @details     The nth root of a complex number is defined as 
 *              \f{eqnarray*}{
 *                  \underline{z} &=& \rho\cdot e^\theta\\
 *                  \sqrt[n]{\underline{z}} &=& \rho^{\frac{1}{n}}\cdot e^{\theta / n}
 *              \f}
 *              where \f$\rho\f$ denotes the absolute value and 
 *              \f$\theta\f$ the argument of the complex number \f$\underline{z}\f$
 *
 * @param[in]   c A complex number, that can represent a real number by specifing 
 *              only the real part of the complex number.
 * @param[in]   n The \f$n\f$-th root is calculated.
 * @tparam      eT An element type which represents a number that provides all common
 *              mathmatical operations.
 *
 * @return      The \f$n\f$-th root as an complex numbers.
 *
 * @since       0.0.1
 *
 * @author      Denis-Michael Lux <denis.lux@icloud.com>
 * @date        12.01.15
 *
 * @ingroup     complex
 */
template< typename eT >
inline
complex< eT > cx_ntrt(const complex< eT >& c, const int n)
{
    complex< eT > com;
    
    eT r = c.abs();
    eT t = c.arg();
    
    if (is_double< eT >::value == true)
    {
        com.polar(pow( r, 1/static_cast< double >(n)), t/static_cast< double >(n));
    }
    else if (is_float< eT >::value == true)
    {
        com.polar(powf(r, 1/static_cast< float >(n)), t/static_cast< float >(n));
    }
    else if (is_ldouble< eT >::value == true)
    {
        com.polar(powl(r, 1/static_cast< long double >(n)), t/static_cast< long double >(n));
    }
    else
    {
        com.polar(static_cast< eT >(pow(r, 1/static_cast< double >(n))), static_cast< eT >(t/static_cast< double >(n)));
    }
    
    return com;
}

#endif
