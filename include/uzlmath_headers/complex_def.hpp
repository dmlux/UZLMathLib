//
//  comlex_def.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 12.01.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_comlex_def_hpp
#define uzlmath_comlex_def_hpp

/*!
 * @brief           Destructor for a complex number.
 * @details         Frees allocated memory
 */
template<typename eT>
inline
complex< eT >::~complex()
{}

/*!
 * @brief           Default constructor for creating a complex number.
 * @details         A default complex number contains only a real value
 *                  of 0 and a imaginary value of 0.
 */
template<typename eT>
inline
complex< eT >::complex()
    : re(0)
    , im(0)
{}



/*!
 * @brief           Constructs a complex number with the same real and imaginary
 *                  value.
 * @details         Constructs a complex number which has the same value for the
 *                  real and the imaginary value.
 * 
 * @param[in]       real_imag The real and complex value.
 *
 */
template<typename eT>
inline
complex< eT >::complex(const eT real_imag)
    : re(real_imag)
    , im(real_imag)
{}



/*!
 * @brief           Constructs a complex number with the given value for the real
 *                  value and the given value for the imaginary value.
 * @details         Constructs a complex number that has a real value of the value
 *                  that is given by the parameter _real_ and a imaginary value 
 *                  given by the parameter _imag_
 *
 * @param[in]       real The real value of the complex number that is supposed to 
 *                  be constructed.
 * @param[in]       imag The imaginary value of the complex number that is supposed
 *                  to be constructed.
 *
 */
template<typename eT>
inline
complex< eT >::complex(const eT real, const eT imag)
    : re(real)
    , im(imag)
{}


/*!
 * @brief           The copy constructor for a complex number.
 * @details         Constructs a complex number that is a copy of the given complex
 *                  number.
 * 
 * @param[in]       c The complex number that is supposed to be copied.
 */
template<typename eT>
inline
complex< eT >::complex(const complex< eT >& c)
    : re(c.re)
    , im(c.im)
{}

/*!
 * @brief           Returns the absolute value of the current complex number.
 * @details         The absolute value of a complex number \f$\underline{z}\f$ 
 *                  is defined as the euclidian norm of real and imaginary.
 *                  \f[
 *                      |\underline{z}| = \sqrt{\Re{(z)}^2 + \Im{(z)}^2}
 *                  \f]
 *
 * @return          The absolute value of the current complex number as defined
 *                  in the detailed description.
 */
template<typename eT>
inline
const eT complex< eT >::abs() const
{
    if (is_double< eT >::value == true)
    {
        return sqrt(re * re + im * im);
    }
    else if (is_float< eT >::value == true)
    {
        return sqrtf(re * re + im * im);
    }
    else if (is_ldouble< eT >::value == true)
    {
        return sqrtl(re * re + im * im);
    }
    else
    {
        return sqrtf(static_cast< float >(re * re) + static_cast< float >(im * im));
    }
}

/*!
 * @brief           Returns the argument of the current complex number.
 * @details         The argument of a complex number \f$\underline{z} = x+yi\f$
 *                  is defined as \f$atan2(z)\f$
 *                  \f[
 *                      \varphi=\arg(\underline{z})=\begin{cases}
 *                          \arctan(\frac{y}{x}) & \mathrm{if}\;x>0\\
 *                          \arctan(\frac{y}{x}) + \pi & \mathrm{if }\;x<0\;\mathrm{ and }\;y\geq 0\\
 *                          \arctan(\frac{y}{x}) - \pi & \mathrm{if }\;x<0\;\mathrm{ and }\;y<0\\
 *                          \frac{\pi}{2} & \mathrm{if }\;x=0\;\mathrm{ and }\;y>0\\
 *                          -\frac{\pi}{2} & \mathrm{if }\;x=0\;\mathrm{ and }\;y<0\\
 *                          \mathrm{indeterminate} & \mathrm{if }\;x=0\;\mathrm{ and }\;y=0
 *                      \end{cases}
 *                  \f]
 *
 * @return          The argument of the current complex number
 */
template<typename eT>
inline
const eT complex< eT >::arg() const
{
    if (is_double< eT >::value == true)
    {
        return atan2(im, re);
    }
    else if (is_float< eT >::value == true)
    {
        return atan2f(im, re);
    }
    else if (is_ldouble< eT >::value == true)
    {
        return atan2l(im, re);
    }
    else
    {
        return atan2f(static_cast< float >(im), static_cast< float >(re));
    }
}

/*!
 * @brief           Returns the norm of the current complex number.
 * @details         The norm of a complex number \f$\underline{z}\f$ is
 *                  defined as
 *                  \f[
 *                      \Re(\underline{z})^2 + \Im(\underline{z})^2
 *                  \f]
 *
 * @return          The norm of the current complex number.
 */
template<typename eT>
inline
const eT complex< eT >::norm() const
{
    return re * re + im * im;
}



/*!
 * @brief           Returns the complexe conjugate of the current complex number.
 * @details         The complex conjugate of a complex number \f$\underline{z}\f$
 *                  is defined as
 *                  \f[
 *                      \rho e^{i\theta\;*} = \rho e^{-i\theta}
 *                  \f]
 */
template<typename eT>
inline
void complex< eT >::conj()
{
    im *= -1;
}

/*!
 * @brief           Constructs a complex number by polar coordinates.
 * @details         The complex number \f$\underline{z}\f$ which is defined by polar 
 *                  coordinates can be expressed in real and imaginary parts by applying
 *                  \f{eqnarray*}{
 *                      \Re(\underline{z}) &=& \rho\cdot \cos(\theta)\\
 *                      \Im(\underline{z}) &=& \rho\cdot \sin(\theta)
 *                  \f}
 *                  where \f$\rho\f$ is the absolute value of the complex number in polar
 *                  coordinates and \f$\theta\f$ is the argument of the copmlex number
 *                  as a radiant angle.
 *
 * @param[in]       rho The absolute value of the complex number.
 * @param[in]       theta The argument of the complex number as radiant angle.
 */
template<typename eT>
inline
void complex< eT >::polar(const eT& rho, const eT& theta)
{
    if (is_double< eT >::value == true)
    {
        re = rho * cos(theta);
        im = rho * sin(theta);
    }
    else if (is_float< eT >::value == true)
    {
        re = rho * cosf(theta);
        im = rho * sinf(theta);
    }
    else if (is_ldouble< eT >::value == true)
    {
        re = rho * cosl(theta);
        im = rho * sinl(theta);
    }
    else
    {
        re = static_cast< eT >(rho * cosf(static_cast< float >(theta)));
        im = static_cast< eT >(rho * cosf(static_cast< float >(theta)));
    }
}

/*!
 * @brief           Addition operator for two complex numbers.
 * @details         Adds two complex nubers by adding the real values
 *                  and the imaginary values of both complex numbers.
 *
 * @param[in]       rhs Complex number on the right handside of the addition
 *                  operator
 * 
 * @return          The result of the addition.
 */
template<typename eT>
inline
complex< eT > complex< eT >::operator+(const complex< eT >& rhs)
{
    return complex< eT >(re + rhs.re, im + rhs.im);
}

/*!
 * @brief           Subtraction operator for two complex numbers.
 * @details         Subtracts the given complex number from the current
 *                  complex number by subtracting the real and imaginary values 
 *                  of the current complex number by the real and imaginary 
 *                  values of the given complex number.
 * 
 * @param[in]       rhs The complex number that is supposed to be subtracted from
 *                  the current complex number.
 *
 * @return          A new complex number containing the result.
 */
template<typename eT>
inline
complex< eT > complex< eT >::operator-(const complex< eT >& rhs)
{
    return complex< eT >(re - rhs.re, im - rhs.im);
}

/*!
 * @brief           Multiplication operator for two complex numbers.
 * @details         Multiplies the given complex number to the current complex
 *                  number by applying the following calculation
 *                  \f{eqnarray*}{
 *                      \underline{z} &=& \underline{x}\cdot \underline{y}\\
 *                      \underline{z} &=& (\Re(\underline{x})\cdot\Re(\underline{y}) 
 *                      - \Im(\underline{x})\cdot\Im(\underline{y}))+(\Im(\underline{x})\cdot\Re(\underline{y}) 
 *                      + \Re(\underline{x})\cdot\Im(\underline{y}))i
 *                  \f}
 *                  where \f$\underline{x},\;\underline{y},\;\underline{z}\f$
 *                  denoting complex numbers.
 *
 * @param[in]       rhs The complex number that is supposed to be multiplied to
 *                  the current complex number.
 *
 * @return          A new complex number containing the multiplication result.
 */
template<typename eT>
inline
complex< eT > complex< eT >::operator*(const complex< eT >& rhs)
{
    return complex< eT >(re * rhs.re - im * rhs.im,   // real value of this number
                         im * rhs.re + re * rhs.im);  // imag value of this number
}

/*!
 * @brief           Division operator for two complex numbers.
 * @details         Divides the given complex number from the current complex
 *                  number by applying the following calculation
 *                  \f{eqnarray*}{
 *                      \underline{z} &=& \underline{x}\cdot \underline{y}\\
 *                      \underline{z} &=& \left(\frac{\Re(\underline{x})\cdot\Re(\underline{y}) 
 *                          - \Im(\underline{x})\cdot\Im(\underline{y})}{\Re(\underline{y})^2 
 *                          + \Im(\underline{y})^2}\right) + \left(\frac{\Im(\underline{x})\cdot\Re(\underline{y})
 *                          +\Re(\underline{x})\cdot\Im(\underline{y})}{\Re(\underline{y})^2 
 *                          + \Im(\underline{y})^2}\right)i
 *                  \f}
 *                  where \f$\underline{x},\;\underline{y},\;\underline{z}\f$
 *                  denoting complex numbers.
 *
 * @param[in]       rhs The complex number that is supposed to be devided by.
 *
 * @return          A new complex number containing the division result.
 */
template<typename eT>
inline
complex< eT > complex< eT >::operator/(const complex< eT >& rhs)
{
    return complex< eT >((re * rhs.re + im * rhs.im) / rhs.norm(),    // real value
                         (im * rhs.re - re * rhs.im) / rhs.norm());   // imag value
}

/*!
 * @brief           Power operator for two complex numbers.
 * @details         Calculates the \f$n\f$-th power of the current complex number
 *                  where \f$n\f$ is the given exponent on the right handside of
 *                  the power operator. The power of the complex number is calculated
 *                  as follows.
 *                  \f{eqnarray*}{
 *                      \underline{z}^n = \rho^n\cdot e^{\theta\cdot n}
 *                  \f}
 *                  where \f$\underline{z}\f$ is a complex number.
 *
 * @param[in]       rhs The exponent of the power operator.
 * 
 * @return          A new complex number containing the result of the
 *                  power operator.
 */
template<typename eT>
inline
complex< eT > complex< eT >::operator^(const int rhs)
{
    complex< eT > c;
    eT t = arg() * rhs;
    eT r;
    
    if (is_double< eT >::value == true)
    {
        r = pow( abs(), static_cast< double >(rhs));
    }
    else if (is_float< eT >::value == true)
    {
        r = powf(abs(), static_cast< float >(rhs));
    }
    else if (is_ldouble< eT >::value == true)
    {
        r = powl(abs(), static_cast< long double >(rhs));
    }
    else
    {
        r = static_cast< eT >( pow(static_cast< double >(abs()), static_cast< double >(rhs)) );
    }
    
    c.polar(r, t);
    return c;
}

/*!
 * @brief           Addition operator for an integer value and a complex number.
 * @details         Calculates the sum of an integer value and a complex number
 *                  by casting the int to a complex number and performing the
 *                  addition of two complex numbers.
 *
 * @param[in]       lhs The integer value on the left handside of the addition
 *                  operator.
 * @param[in]       rhs The complex value on the right handside of the addition
 *                  operator.
 *
 * @return          A new complex number containing the result of the addition.
 * 
 * @sa              complex::operator+
 */
template<typename eT>
complex< eT > operator+(int lhs, complex< eT > rhs)
{
    return complex< eT >(lhs, 0) + rhs;
}

/*!
 * @brief           Addition operator for an non-complex value and a complex number.
 * @details         Calculates the sum of an non-complex value and a complex number
 *                  by casting the non-complex to a complex number and performing the
 *                  addition of two complex numbers.
 *
 * @param[in]       rhs The non-complex value on the right handside of the addition
 *                  operator.
 *
 * @return          A new complex number containing the result of the addition.
 *
 * @sa              complex::operator+
 */
template<typename eT>
template<typename T>
inline
complex< eT > complex< eT >::operator+(const T& rhs)
{
    complex< eT > c(rhs, 0);
    return *this + c;
}

/*!
 * @brief           Subtraction operator for an non-complex value and a complex number.
 * @details         Calculates the difference of an non-complex value and a complex number
 *                  by casting the non-complex to a complex number and performing the
 *                  subtraction of two complex numbers.
 *
 * @param[in]       rhs The non-complex value on the right handside of the subtraction
 *                  operator.
 *
 * @return          A new complex number containing the result of the subtraction.
 *
 * @sa              complex::operator-
 */
template<typename eT>
template<typename T>
inline
complex< eT > complex< eT >::operator-(const T& rhs)
{
    complex< eT > c(rhs, 0);
    return *this - c;
}

/*!
 * @brief           Multiplication operator for an non-complex value and a complex number.
 * @details         Calculates the product of an non-complex value and a complex number
 *                  by casting the non-complex to a complex number and performing the
 *                  multiplication of two complex numbers.
 *
 * @param[in]       rhs The non-complex value on the right handside of the multiplication
 *                  operator.
 *
 * @return          A new complex number containing the result of the multiplication.
 *
 * @sa              complex::operator*
 */
template<typename eT>
template<typename T>
inline
complex< eT > complex< eT >::operator*(const T& rhs)
{
    complex< eT > c(rhs, 0);
    return *this * c;
}

/*!
 * @brief           Division operator for an non-complex value and a complex number.
 * @details         Calculates the quotient of an non-complex value and a complex number
 *                  by casting the non-complex to a complex number and performing the
 *                  division of two complex numbers.
 *
 * @param[in]       rhs The non-complex value on the right handside of the division
 *                  operator.
 *
 * @return          A new complex number containing the result of the division.
 *
 * @sa              complex::operator/
 */
template<typename eT>
template<typename T>
inline
complex< eT > complex< eT >::operator/(const T& rhs)
{
    complex< eT > c(rhs, 0);
    return *this / c;
}



/*!
 * @brief           The addition assignment operator for two complex numbers.
 * @details         Adds a complex number to the current complex number and 
 *                  storing the result in the current complex number by 
 *                  overwriting the contents of the current complex number.
 *
 * @param[in]       rhs The complex number on the right handside of the 
 *                  addition assignment operator.
 * 
 * @return          The reference to the current complex number
 *
 * @sa              For further information on how to add two complex numbers
 *                  see complex::operator+
 */
template<typename eT>
inline
const complex< eT >& complex< eT >::operator+=(const complex< eT >& rhs)
{
    re += rhs.re;
    im += rhs.im;
    return *this;
}

/*!
 * @brief           The subtraction assignment operator for two complex numbers.
 * @details         Subtracts a complex number from the current complex number and
 *                  storing the result in the current complex number by
 *                  overwriting the contents of the current complex number.
 *
 * @param[in]       rhs The complex number on the right handside of the
 *                  subtraction assignment operator.
 *
 * @return          The reference to the current complex number
 *
 * @sa              For further information on how to subtract two complex numbers
 *                  see complex::operator-
 */
template<typename eT>
inline
const complex< eT >& complex< eT >::operator-=(const complex< eT >& rhs)
{
    re -= rhs.re;
    im -= rhs.im;
    return *this;
}

/*!
 * @brief           The multiplication assignment operator for two complex numbers.
 * @details         Multiplies a complex number to the current complex number and
 *                  storing the result in the current complex number by
 *                  overwriting the contents of the current complex number.
 *
 * @param[in]       rhs The complex number on the right handside of the
 *                  multiplication assignment operator.
 *
 * @return          The reference to the current complex number
 *
 * @sa              For further information on how to multiply two complex numbers
 *                  see complex::operator*
 */
template<typename eT>
inline
const complex< eT >& complex< eT >::operator*=(const complex< eT >& rhs)
{
    complex< eT > c = *this * rhs;
    re = c.re;
    im = c.im;
    return *this;
}

/*!
 * @brief           The division assignment operator for two complex numbers.
 * @details         Divides a complex number by a complex number and
 *                  storing the result in the current complex number by
 *                  overwriting the contents of the current complex number.
 *
 * @param[in]       rhs The complex number on the right handside of the
 *                  division assignment operator.
 *
 * @return          The reference to the current complex number
 *
 * @sa              For further information on how to divide two complex numbers
 *                  see complex::operator/
 */
template<typename eT>
inline
const complex< eT >& complex< eT >::operator/=(const complex< eT >& rhs)
{
    complex< eT > c = *this / rhs;
    re = c.re;
    im = c.im;
    return *this;
}

/*!
 * @brief           The power assignment operator for a complex numbers.
 * @details         Calculates the power of a complex number and storing the 
 *                  result in the current complex number by overwriting the 
 *                  contents of the current complex number.
 *
 * @param[in]       rhs Exponent of the power assignment operator.
 *
 * @return          The reference to the current complex number
 *
 * @sa              For further information on how to multiply two complex numbers
 *                  see complex::operator*
 */
template<typename eT>
inline
const complex< eT >& complex< eT >::operator^=(const int rhs)
{
    eT t = arg() * rhs;
    eT r;
    
    if (is_double< eT >::value == true)
    {
        r = pow( abs(), static_cast< double >(rhs));
    }
    else if (is_float< eT >::value == true)
    {
        r = powf(abs(), static_cast< float >(rhs));
    }
    else if (is_ldouble< eT >::value == true)
    {
        r = powl(abs(), static_cast< long double >(rhs));
    }
    else
    {
        r = static_cast< eT >(powf(static_cast< float >(abs()), rhs));
    }
    
    polar(r, t);
    return *this;
}

/*!
 * @brief           The addition assignment operator for a complex number and
 *                  a non-complex number.
 * @details         Adds a non-complex number to the current complex number and
 *                  storing the result in the current complex number by
 *                  overwriting the contents of the current complex number. The
 *                  non-complex number gets casted into a complex number and
 *                  is then added to the current complex number.
 *
 * @param[in]       rhs The complex number on the right handside of the
 *                  addition assignment operator.
 *
 * @return          The reference to the current complex number
 *
 * @sa              For further information on how to add two complex numbers
 *                  see complex::operator+
 */
template<typename eT>
template<typename T>
inline
const complex< eT >& complex< eT >::operator+=(const T& rhs)
{
    if (is_complex< T >::value == true)
    {
        *this += rhs;
    }
    else
    {
        complex< eT > c(rhs, 0);
        *this += c;
    }
    return *this;
}

/*!
 * @brief           The subtraction assignment operator for a complex number and
 *                  a non-complex number.
 * @details         Subtracts a non-complex number from the current complex number and
 *                  storing the result in the current complex number by
 *                  overwriting the contents of the current complex number. The
 *                  non-complex number gets casted into a complex number and
 *                  is then subtracted from the current complex number.
 *
 * @param[in]       rhs The complex number on the right handside of the
 *                  subtraction assignment operator.
 *
 * @return          The reference to the current complex number
 *
 * @sa              For further information on how to subtract two complex numbers
 *                  see complex::operator-
 */
template<typename eT>
template<typename T>
inline
const complex< eT >& complex< eT >::operator-=(const T& rhs)
{
    if (is_complex< T >::value == true)
    {
        *this -= rhs;
    }
    else
    {
        complex< eT > c(rhs, 0);
        *this -= c;
    }
    return *this;
}

/*!
 * @brief           The mutliplication assignment operator for a complex number and
 *                  a non-complex number.
 * @details         Multiplies a non-complex number from the current complex number and
 *                  storing the result in the current complex number by
 *                  overwriting the contents of the current complex number. The
 *                  non-complex number gets casted into a complex number and
 *                  is then multiplied with the current complex number.
 *
 * @param[in]       rhs The complex number on the right handside of the
 *                  mutliplication assignment operator.
 *
 * @return          The reference to the current complex number
 *
 * @sa              For further information on how to multiply two complex numbers
 *                  see complex::operator*
 */
template<typename eT>
template<typename T>
inline
const complex< eT >& complex< eT >::operator*=(const T& rhs)
{
    if (is_complex< T >::value == true)
    {
        *this *= rhs;
    }
    else
    {
        complex< eT > c(rhs, 0);
        *this *= c;
    }
    return *this;
}

/*!
 * @brief           The division assignment operator for a complex number and
 *                  a non-complex number.
 * @details         Divides a non-complex number from the current complex number and
 *                  storing the result in the current complex number by
 *                  overwriting the contents of the current complex number. The
 *                  non-complex number gets casted into a complex number and
 *                  then the division is performed.
 *
 * @param[in]       rhs The complex number on the right handside of the
 *                  division assignment operator.
 *
 * @return          The reference to the current complex number
 *
 * @sa              For further information on how to divide two complex numbers
 *                  see complex::operator/
 */
template<typename eT>
template<typename T>
inline
const complex< eT >& complex< eT >::operator/=(const T& rhs)
{
    if (is_complex< T >::value == true)
    {
        *this /= rhs;
    }
    else
    {
        complex< eT > c(rhs, 0);
        *this /= c;
    }
    return *this;
}


/*- Additional operators -*/


/*!
 * @brief           Multiplication operator for an non-complex value and a complex number.
 * @details         Calculates the product of an non-complex value and a complex number
 *                  by casting the non-complex to a complex number and performing the
 *                  multiplication of two complex numbers.
 *
 * @param[in]       rhs The non-complex value on the right handside of the multiplication
 *                  operator.
 *
 * @return          A new complex number containing the result of the multiplication.
 *
 * @sa              complex::operator*
 */
template<typename eT, typename T>
complex< eT > operator*(const T& lhs, complex< eT > rhs)
{
    return complex< eT >(lhs, 0) * rhs;
}

/*!
 * @brief           Division operator for an non-complex value and a complex number.
 * @details         Calculates the quotient of an non-complex value and a complex number
 *                  by casting the non-complex to a complex number and performing the
 *                  division of two complex numbers.
 *
 * @param[in]       rhs The non-complex value on the right handside of the division
 *                  operator.
 *
 * @return          A new complex number containing the result of the division.
 *
 * @sa              complex::operator/
 */
template<typename eT, typename T>
complex< eT > operator/(const T& lhs, complex< eT > rhs)
{
    return complex< eT >(lhs, 0) / rhs;
}

/*!
 * @brief           Addition operator for an non-complex value and a complex number.
 * @details         Calculates the sum of an non-complex value and a complex number
 *                  by casting the non-complex to a complex number and performing the
 *                  addition of two complex numbers.
 *
 * @param[in]       rhs The non-complex value on the right handside of the addition
 *                  operator.
 *
 * @return          A new complex number containing the result of the addition.
 *
 * @sa              complex::operator+
 */
template<typename eT, typename T>
complex< eT > operator+(const T& lhs, complex< eT > rhs)
{
    return complex< eT >(lhs, 0) + rhs;
}

/*!
 * @brief           Subtraction operator for an non-complex value and a complex number.
 * @details         Calculates the difference of an non-complex value and a complex number
 *                  by casting the non-complex to a complex number and performing the
 *                  subtraction of two complex numbers.
 *
 * @param[in]       rhs The non-complex value on the right handside of the subtraction
 *                  operator.
 *
 * @return          A new complex number containing the result of the subtraction.
 *
 * @sa              complex::operator-
 */
template<typename eT, typename T>
complex< eT > operator-(const T& lhs, complex< eT > rhs)
{
    return complex< eT >(lhs, 0) - rhs;
}

/*!
 * @brief           Outstream operator overload to print the complex number in a nice
 *                  format to the outstream.
 * @details         The overload for the outstream operator can be used like this
 *                  @code
 *                      complex<double> z(2.5, 3.2);
 *                      std::cout << "z = " << z << std::endl;
 *                      
 *                      // The output looks as follows:
 *                      // z = +2.5000e+00+3.2000e+00i
 *                  @endcode
 *                  The number is printed in scientific format.
 *
 * @param[in]       o The outstream object.
 * @param[in]       c The complex number that is supposed to be printed.
 * @tparam          eT The template type of the given complex number.
 *
 * @return          The reference to the outstream object.
 */
template<typename eT>
std::ostream& operator<<(std::ostream& o, const complex< eT >& c)
{
    std::ios::fmtflags f( std::cout.flags() );
    o << std::scientific << std::setprecision(4);
    
    o << std::noshowpos;
    o << c.re;
    o << (c.im < 0 ? " - " : " + ") << (c.im == 0 ?  0 : UZL_ABS(c.im)) << "i";
    
    std::cout.flags( f );
    return o << std::fixed;
}

#endif
