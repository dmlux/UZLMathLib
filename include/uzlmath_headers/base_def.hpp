//
//  base_def.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 21.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_base_def_hpp
#define uzlmath_base_def_hpp

template< typename derived >
inline
const derived& base< derived >::get_ref() const
{
    return static_cast< const derived& >(*this);
}

template< typename T1, typename T2 >
inline
glue< T1, T2 > operator+(const base< T1 >& A, const base< T2 >& B)
{
    return glue< T1, T2 >(A.get_ref(), B.get_ref());
}

#endif
