//
//  random.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 04.10.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_random_hpp
#define UZLMathLib_random_hpp

UZLMATH_BEGIN

/*!
 * @brief       A context for random number creation
 */
template< typename pod_type, typename derived >
struct randctx
{
    inline const derived& get_ref() const;
};

template< typename pod_type, typename derived >
inline
const derived& randctx< pod_type, derived >::get_ref() const
{
    return static_cast< const derived& >(*this);
}

UZLMATH_END

#endif /* random.hpp */
