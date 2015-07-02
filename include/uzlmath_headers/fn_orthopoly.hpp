//
//  fn_orthopoly.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 12.01.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_fn_orthopoly_hpp
#define uzlmath_fn_orthopoly_hpp

/*!
 * @brief   The polynomials namespace contains functions that defining
 *          polynomials and functions that are used in analysis.
 *
 * @since   0.0.1
 *
 * @author  Denis-Michael Lux <denis.lux@icloud.com>
 * @date    12.01.15
 */
namespace orthopoly
{
    // Legendre polynome for x of degree n
    auto legendre(const int& n, const double& x) -> const double;
}

#endif
