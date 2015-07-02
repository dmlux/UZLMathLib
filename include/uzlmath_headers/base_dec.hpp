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

/*!
 * @brief           Collection of base class classes and functions.
 * @defgroup        base Base class
 * @{
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

#endif
