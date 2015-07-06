//
//  factorial.cpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 12.01.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_factorial_cpp
#define uzlmath_factorial_cpp

#include <uzlmath>

UZLMATH_BEGIN

factorial::factorial()
    : value(0)
    , exponents(nullptr)

{}

factorial::factorial(const int& number)
    : value(number)
{
    // allocate dynamic memory
    exponents = new int[number];
    
    // set values in exponents
    unsigned int i;
    for (i = 0; i < number; ++i)
    {
        exponents[i] = 1;
    }
}

factorial::factorial(const factorial& factorial)
    : value(factorial.value)
{
    exponents = new int[factorial.value];
    memcpy(exponents, factorial.exponents, sizeof(int) * factorial.value);
}

factorial::factorial(factorial&& factorial)
    : value(factorial.value)
{
    int* tmp            = exponents;
    exponents           = factorial.exponents;
    factorial.exponents = tmp;
}

factorial::~factorial()
{
    delete [] exponents;
}

const
factorial& factorial::operator=(const factorial& rhs)
{
    if ( this == &rhs )
    {
        return *this;
    }
    
    delete [] exponents;
    
    exponents = new int[rhs.value];
    memcpy(exponents, rhs.exponents, sizeof(int) * rhs.value);
    
    return *this;
}

const
factorial& factorial::operator=(factorial&& rhs)
{
    if ( this == &rhs )
    {
        return *this;
    }
    
    delete [] exponents;
    
    exponents = rhs.exponents;
    rhs.exponents = nullptr;
    
    return *this;
}

factorial factorial::operator*(const factorial& rhs)
{
    // create result
    int n_bigger  = value < rhs.value  ? rhs.value  : value;
    int n_smaller = value < rhs.value  ? value      : rhs.value;
    
    int* bigger   = n_bigger  == value ? exponents  : rhs.exponents;
    int* smaller  = n_smaller == value ? exponents  : rhs.exponents;
    
    factorial result(n_bigger);
    memcpy(result.exponents, bigger, sizeof(int) * n_bigger);
    
    unsigned int i;
    for (i = 0; i < n_smaller; ++i)
    {
        result.exponents[i] += smaller[i];
    }
    
    return result;
}

factorial factorial::operator/(const factorial& rhs)
{
    // create result
    int n_bigger = value < rhs.value ? rhs.value : value;
    
    factorial result(n_bigger);
    
    memset(result.exponents, 0, sizeof(int) * n_bigger);
    memcpy(result.exponents, exponents, sizeof(int) * value);
    
    unsigned int i;
    for (i = 0; i < rhs.value; ++i)
    {
        result.exponents[i] -= rhs.exponents[i];
    }
    
    return result;
}

double factorial::eval() const
{
    double result = 1;
    
    unsigned int i;
    for (i = 0; i < value; ++i)
    {
        if (exponents[i] == 1)
        {
            result *= i + 1;
        }
        else if (exponents[i] != 0)
        {
            result *= pow(i + 1, exponents[i]);
        }
    }
    
    return result;
}

UZLMATH_END
    
#endif
