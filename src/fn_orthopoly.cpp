//
//  fn_orthopoly.cpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 12.01.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_fn_orthopoly_cpp
#define UZLMathLib_fn_orthopoly_cpp

#include <uzlmath>

/*!
 * @brief   The polynomials namespace contains functions that defining
 *          polynomials and functions that are used in analysis.
 *
 * @since   0.0.1
 *
 * @author  Denis-Michael Lux <denis.lux@icloud.com>
 * @date    12.01.15
 */
UZLMATH_NAMESPACE(orthoPoly)

/*!
 * @brief       Computes the value of the Legendre polynomial \f$P_n(x) = \frac{(2n-1)x}{n}P_{n-1}(x)-(1-\frac{1}{n})P_{n-2}(x)\f$
 * @details     The legendre functions are solutions to Legendre's differential equation
 *              \f{eqnarray*}{
 *                  \frac{\mathrm d}{\mathrm d x}\left[(1-x^2)\;\frac{\mathrm d}{\mathrm d x}P_n(x)\right] + n(n+1)P_n(x) = 0
 *              \f}
 *              The implementation of this function is based on the three-term recurrence
 *              that is shown above. The recurrence is solved by use of dynamic programming
 *              to provide an \f$\mathcal{O}(n)\f$ algorithm to calculate the result.
 *
 *              One important property of the Legendre polynomials is that they are
 *              orthogonal with respect to the \f$L^2\f$ inner product on the interval
 *              \f$-1\leq x\leq 1\f$. Additionally they are symmetric and antisymmetric
 *
 * @param[in]   n Degree \f$n\f$ of the polynome
 * @param[in]   x Value \f$x\f$ of the polynome
 *
 * @author      Denis-Michael Lux <denis.lux@icloud.com>
 * @date        12.01.15
 *
 * @since       0.0.1
 */
double legendre(const int& n, const double& x)
{
    // illegal values
    if (n < 0)
    {
        uzlmath_warning("%s", "degree of legendre polynomial is negative!");
        return 0.0;
    }
    
    // intialize three-term recurrence
    double preprev  = 1;
    double prev     = x;
    
    // start with 2 because first two values are predefined with
    // 1 and x
    for (int i = 2; i <= n; ++i)
    {
        // make copy of the previous value
        double tmp  = prev;
        
        // shift values
        prev        = (2. * i - 1.) / i * x * prev - (i - 1.) / i * preprev;
        preprev     = tmp;
    }
    
    // return previous value which represents value of the legendre polynomial
    return prev;
}

/*!
 * @brief       Calculates the associated legendre polynomials needed
 *              To calculate the spherical harmonics
 */
double assoc_legendre(const int& l, const int& m, const double& x)
{
    if ( l < m || l < abs(m) )
    {
        uzlmath_warning("%s", "Order l of associated legendre polynomial is greater degree |m|!");
        return 0.0;
    }
    
    // indices
    int i, pos_m = abs(m);
    
    // get signs of polyonmial
    double sign = 1;
    
    // correct sign for negative m. If m is negative
    // The aboslute value of m is taken for the poly-
    // nomial calculation and a scalar is muliplied
    // with the polyonmial. The factor is given by
    //
    //  (-1)^m * (l - m)!/(l + m)!
    //
    if ( m < 0 )
    {
        // calculate factorials
        double lmm = 1;
        for (i = 1; i <= l - pos_m; ++i)
        {
            lmm *= i;
        }
        
        double lpm = 1;
        for (i = 1; i <= l + pos_m; ++i)
        {
            lpm *= i;
        }
        
        // update sign variable
        sign *= (pos_m & 1 ? -1.0 : 1.0) * lmm / lpm;
    }
    
    // *** Calculate PREPREV ***
    //
    // first get double factorial of (2m - 1). Then
    // calculate the last scaling factor which is
    // given by (1 - x^2)^(l / 2). At last multiply
    // the scaling factor on the preprev variable
    // which contains the product of the other terms.
    // The initial PREPREV is defined as
    //
    //  (-1)^l * (2l - 1)!!(1 - x^2) ^ (l / 2)
    //
    double preprev = pos_m & 1 ? -1.0 : 1.0, scale = 1;
    
    for (i = 1; i <= 2 * pos_m - 1; i += 2)
    {
        preprev *= i;
    }
    
    for (i = 0; i < pos_m; ++i)
    {
        scale *= 1.0 - x * x;
    }
    
    preprev *= sqrt( scale );
    
    // *** Calculate PREV ***
    //
    // Calculate the previous recurrence value
    // which is defined as
    //
    //  PREV = x * (2l + 1) * PREPREV
    //
    double prev = x * (2.0 * pos_m + 1.0) * preprev;
    
    
    // *** Return values ***
    //
    // if pos_m == l we already have calculated the value in preprev
    if ( pos_m == l )
    {
        return sign * preprev;
    }
    
    // if pos_m == l + 1 we already have caculated the value in prev
    else if ( pos_m == l + 1 )
    {
        return sign * prev;
    }
    
    // *** Calculate three-term-recurrence ***
    //
    // now iterate until from m+2 until l is reaced
    else
    {
        // start three-term-recurrence
        for (i = pos_m + 2; i <= l; ++i)
        {
            // store current previous value
            double tmp = prev;
            
            // shift values
            prev = x * (2.0 * i - 1.0) / (i - pos_m) * prev - (i + pos_m - 1.0) / (i - pos_m) * preprev;
            preprev = tmp;
        }
        
        // return recurrence value
        return sign * prev;
    }
}

UZLMATH_NAMESPACE_END

#endif
