//
//  fn_cx_sqrt.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 12.01.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_fn_cx_sqrt_hpp
#define UZLMathLib_fn_cx_sqrt_hpp

UZLMATH_BEGIN

/*!
 * @brief       Calculates the square root of a given real number. And
 *              allows to calculate the square root of negative real numbers.
 *
 * @param[in]   n A real number that can be positive as well as negative.
 * @tparam      eT An element type which represents a number that provides all common
 *              mathmatical operations.
 * @return      The square root as an complex numbers. To provide the
 *              possibility of calculating the square root of negative
 *              values the result has to be complex. If the value is
 *              positive the resulting complex number has only a real
 *              part.
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
complex< eT > cx_sqrt(const eT n)
{
    complex< eT > c;
    
    if (n < 0)
    {
        if (is_double< eT >::value == true)
        {
            c.im = sqrt(-n);
        }
        else if (is_float< eT >::value == true)
        {
            c.im = sqrtf(-n);
        }
        else if (is_ldouble< eT >::value == true)
        {
            c.im = sqrtl(-n);
        }
        else
        {
            c.im = static_cast< eT >(sqrt(static_cast< double >(-n)));
        }
    }
    else
    {
        if (is_double< eT >::value == true)
        {
            c.re = sqrt(n);
        }
        else if (is_float< eT >::value == true)
        {
            c.re = sqrtf(n);
        }
        else if (is_ldouble< eT >::value == true)
        {
            c.re = sqrtl(n);
        }
        else
        {
            c.re = static_cast< eT >(sqrt(static_cast< double >(n)));
        }
    }
    
    return c;
}

/*!
 * @brief       Calculates the square root of a given complex number.
 *
 * @tparam      eT An element type which represents a number that provides all common
 *              mathmatical operations.
 * @param[in]   c A complex number.
 * @return      The square root as an complex numbers.
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
complex< eT > cx_sqrt(const complex< eT >& c)
{
    if (c.im == 0)
    {
        return c_sqrt(c.re);
    }
    
    complex< eT > com;
    short sgn = (c.im < 0 ? -1 : (c.im > 0 ? 1 : 0));
    
    if (is_double< eT >::value == true)
    {
        com.re =       sqrt(( c.re + sqrt(c.re * c.re + c.im * c.im)) / 2.)  ;
        com.im = sgn * sqrt((-c.re + sqrt(c.re * c.re + c.im * c.im)) / 2.)  ;
    }
    else if (is_float< eT >::value == true)
    {
        com.re =       sqrtf(( c.re + sqrtf(c.re * c.re + c.im * c.im)) / 2.);
        com.im = sgn * sqrtf((-c.re + sqrtf(c.re * c.re + c.im * c.im)) / 2.);
    }
    else if (is_ldouble< eT >::value == true)
    {
        com.re =       sqrtl(( c.re + sqrtl(c.re * c.re + c.im * c.im)) / 2.);
        com.im = sgn * sqrtl((-c.re + sqrtl(c.re * c.re + c.im * c.im)) / 2.);
    }
    else
    {
        com.re = static_cast< eT >(sqrt(static_cast< double >(c.re + sqrt(static_cast< double >(c.re * c.re + c.im * c.im))) / 2));
        com.im = static_cast< eT >(sgn * sqrt(static_cast< double >(-c.re + sqrt(c.re * c.re + c.im * c.im)) / 2));
    }
    
    return com;
}

UZLMATH_END

#endif
