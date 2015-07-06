//
//  fn_orthopoly.cpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 12.01.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_fn_orthopoly_cpp
#define uzlmath_fn_orthopoly_cpp

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
auto legendre(const int& n, const double& x) -> const double
{
    // memory for memorization
    double mem[n + 1];
    
    mem[0] = 1;
    mem[1] = x;
    
    for (int i = 2; i <= n; ++i)
    {
        mem[i] = (2.0 * i - 1.0) / i * x * mem[i - 1] - (i - 1.0) / i * mem[i - 2];
    }
    
    return mem[n];
}

UZLMATH_NAMESPACE_END

#endif
