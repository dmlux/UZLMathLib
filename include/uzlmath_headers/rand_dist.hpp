//
//  rand_dist.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 04.10.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_rand_dist_hpp
#define UZLMathLib_rand_dist_hpp

UZLMATH_BEGIN

template< typename T >
struct uniform_real_distribution< T, typename is_true< is_real_type< T >::value >::type > : randctx< T, uniform_real_distribution< T > >
{
public:
    random_engine engine;   //!< Random engine that should be used
    T min;                  //!< Minimum random value
    T max;                  //!< Maximum random value
};

template< typename T >
struct uniform_int_distribution< T, typename is_true< is_integral_type< T >::value >::type > : randctx< T, uniform_int_distribution< T > >
{
public:
    random_engine engine;   //!< Random engine that should be used
    T min;                  //!< Minimum random value
    T max;                  //!< Maximum random value
};

UZLMATH_END

#endif
