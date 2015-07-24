//
//  base_dec.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 23.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_base_dec_hpp
#define uzlmath_base_dec_hpp

UZLMATH_BEGIN

/*!
 * @brief           Collection of base class classes and functions.
 * @defgroup        base Base class
 * @{
 */

/*!
 * @brief           The base class is a thin wrapper around any possible
 *                  object. 
 * @details         The base class can return return a reference to the object
 *                  it wraps.
 *
 * @since           0.0.1
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            24.07.15
 */
template< typename derived >
struct base
{
    inline const derived& get_ref() const;
};

template< typename T1, typename T2 > inline glue< T1, T2 > operator+(const base< T1 >& A, const base< T2 >& B);

/*!
 * @}
 */

UZLMATH_END

#endif
