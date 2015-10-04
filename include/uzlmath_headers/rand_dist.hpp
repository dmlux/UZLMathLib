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

///*!
// * @brief       A uniform distribution for integral types
// */
//template< typename T >
//struct uniform_int_distribution< T, typename is_true< dist_is_integral< T >::value >::type > : random< T, uniform_int_distribution< T > >
//{
//public:
//    random_engine engine;   //!< Random engine that should be used
//    T min;                 //!< Minimum random value
//    T max;                 //!< Maximum random value
//    
//    inline uniform_int_distribution();
//    inline uniform_int_distribution(const uniform_int_distribution< T >& uid);
//};
//
//template< typename T >
//inline
//uniform_int_distribution< T, typename is_true< dist_is_integral< T >::value >::type >::uniform_int_distribution()
//{
//    engine = random_engine::DEFAULT;
//    min = 0;
//    max = 0;
//}
//
//template< typename T >
//inline
//uniform_int_distribution< T, typename is_true< dist_is_integral< T >::value >::type >::uniform_int_distribution(const uniform_int_distribution< T >& uid)
//{
//    engine = uid.engine;
//    min = uid.min;
//    max = uid.max;
//}

UZLMATH_END

#endif
