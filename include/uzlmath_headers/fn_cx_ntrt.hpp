//
//  fn_cx_ntrt.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 12.01.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_fn_cx_ntrt_hpp
#define UZLMathLib_fn_cx_ntrt_hpp

UZLMATH_BEGIN

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
 * @tparam      T An element type which represents a number that provides all common
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
template< typename T >
inline
complex< T > cx_ntrt(const complex< T >& c, const int n)
{
    complex< T > com;
    
    T r = c.abs();
    T t = c.arg();
    
    if ( same_type< T, double >::value )
    {
        com.polar(pow( r, 1/static_cast< double >(n)), t/static_cast< double >(n));
    }
    else if ( same_type< T, float >::value )
    {
        com.polar(powf(r, 1/static_cast< float >(n)), t/static_cast< float >(n));
    }
    else if ( same_type< T, long double >::value )
    {
        com.polar(powl(r, 1/static_cast< long double >(n)), t/static_cast< long double >(n));
    }
    else
    {
        com.polar(static_cast< T >(pow(r, 1/static_cast< double >(n))), static_cast< T >(t/static_cast< double >(n)));
    }
    
    return com;
}

UZLMATH_END

#endif /* fn_cx_ntrt.hpp */
