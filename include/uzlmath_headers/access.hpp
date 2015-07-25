//
//  access.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 25.07.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_access_hpp
#define UZLMathLib_access_hpp

UZLMATH_BEGIN

class access
{
public:
    
    /*!
     * @brief           Removes the constantness of given element
     */
    template< typename T > uzlmath_inline static T&  rw (const T& x)        { return const_cast< T& >(x);  }
    template< typename T > uzlmath_inline static T*& rwp(const T* const& x) { return const_cast< T*& >(x); };
};

UZLMATH_END

#endif
