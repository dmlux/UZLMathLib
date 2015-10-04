//
//  constants.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 29.09.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_constants_hpp
#define UZLMathLib_constants_hpp

UZLMATH_BEGIN

template< typename T >
class constants
{
public:
    static const T pi; //!< The ratio of a circle's circumference to its diameter
    static const T e;  //!< e is the limit of (1 + 1/n)^n for n to infinity
};

template< typename T >
const T constants< T >::pi = static_cast< T >(3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664);

template< typename T >
const T constants< T >::e = static_cast< T >(2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200305992);

UZLMATH_END

#endif
