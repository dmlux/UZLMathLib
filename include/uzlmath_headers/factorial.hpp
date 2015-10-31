//
//  factorial_dec.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 01.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_factorial_dec_hpp
#define UZLMathLib_factorial_dec_hpp

UZLMATH_BEGIN

/*! 
 * @brief   A C++ implementation of a factorial number.
 * @details It can be used to do calculations with big factorial values. 
 *          The number is represented as a multiplication of powers of 
 *          the number.
 *
 * @since   0.0.1
 *
 * @author  Denis-Michael Lux <denis.lux@icloud.com>
 * @date    01.05.15
 *
 * @todo    Implement power operations and all other common operators.
 * @todo    Write additional coversion overload operators.
 */
class factorial
{
    int value;                  //!< The value of factorial that should be computed
    int* exponents;             //!< The exponents that are used to calculate with the factorials
    
public:
    // mark constructors as explicit to prevent constructing a factorial
    // via copy constructor with double value or other conversion type.
             factorial();
    explicit factorial(const int& number);
             factorial(const factorial& factorial);
             factorial(factorial&& factorial);
            ~factorial();
    
    const    factorial& operator=(const factorial& rhs);
    const    factorial& operator=(factorial&& rhs);
    
             factorial  operator*(const factorial& rhs);
             factorial  operator/(const factorial& rhs);
    
    double   eval() const;
};

UZLMATH_END

#endif /* factorial.hpp */
