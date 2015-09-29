//
//  constants.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 29.09.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef constants_hpp
#define constants_hpp

UZLMATH_BEGIN

template< typename eT >
class constants
{
public:
    static const eT pi; //!< The ratio of a circle's circumference to its diameter
    static const eT e;  //!< e is the limit of (1 + 1/n)^n
};

template< typename eT >
const eT constants< eT >::pi = static_cast< eT >(3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664);

template< typename eT >
const eT constants< eT >::e = static_cast< eT >(2.71828182845904523536028747135266249775724709369995957496696762772407663035354759457138217852516642742746639193200305992);

UZLMATH_END

#endif /* constants_h */
