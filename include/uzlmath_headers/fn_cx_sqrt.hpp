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
 * @tparam      T An element type which represents a number that provides all common
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
template< typename T >
inline
complex< T > cx_sqrt(const T n)
{
    complex< T > c;
    
    if (n < 0)
    {
        if ( same_type< T, double >::value )
        {
            c.im = sqrt(-n);
        }
        else if ( same_type< T, float >::value )
        {
            c.im = sqrtf(-n);
        }
        else if ( same_type< T, long double >::value )
        {
            c.im = sqrtl(-n);
        }
        else
        {
            c.im = static_cast< T >(sqrt(static_cast< double >(-n)));
        }
    }
    else
    {
        if ( same_type< T, double >::value )
        {
            c.re = sqrt(n);
        }
        else if ( same_type< T, float >::value )
        {
            c.re = sqrtf(n);
        }
        else if ( same_type< T, long double >::value )
        {
            c.re = sqrtl(n);
        }
        else
        {
            c.re = static_cast< T >(sqrt(static_cast< double >(n)));
        }
    }
    
    return c;
}

/*!
 * @brief       Calculates the square root of a given complex number.
 *
 * @tparam      T An element type which represents a number that provides all common
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
template< typename T >
inline
complex< T > cx_sqrt(const complex< T >& c)
{
    if (c.im == 0)
    {
        return c_sqrt(c.re);
    }
    
    complex< T > com;
    short sgn = (c.im < 0 ? -1 : (c.im > 0 ? 1 : 0));
    
    if ( same_type< T, double >::value )
    {
        com.re =       sqrt(( c.re + sqrt(c.re * c.re + c.im * c.im)) / 2.)  ;
        com.im = sgn * sqrt((-c.re + sqrt(c.re * c.re + c.im * c.im)) / 2.)  ;
    }
    else if ( same_type< T, float >::value )
    {
        com.re =       sqrtf(( c.re + sqrtf(c.re * c.re + c.im * c.im)) / 2.);
        com.im = sgn * sqrtf((-c.re + sqrtf(c.re * c.re + c.im * c.im)) / 2.);
    }
    else if ( same_type< T, long double >::value )
    {
        com.re =       sqrtl(( c.re + sqrtl(c.re * c.re + c.im * c.im)) / 2.);
        com.im = sgn * sqrtl((-c.re + sqrtl(c.re * c.re + c.im * c.im)) / 2.);
    }
    else
    {
        com.re = static_cast< T >(sqrt(static_cast< double >(c.re + sqrt(static_cast< double >(c.re * c.re + c.im * c.im))) / 2));
        com.im = static_cast< T >(sgn * sqrt(static_cast< double >(-c.re + sqrt(c.re * c.re + c.im * c.im)) / 2));
    }
    
    return com;
}

UZLMATH_END

#endif /* fn_cx_sqrt.hpp */
