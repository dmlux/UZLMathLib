//
//  fn_wigner_d.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 12.01.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_fn_wigner_d_hpp
#define uzlmath_fn_wigner_d_hpp

/*!
 * @brief   The **Wigner** namespace contains the Wigner functions and its cousin like the Wigner D-function
 *          and the \f$L^2\f$-normalized Wigner d-function.
 *
 * @since   0.0.1
 *
 * @author  Denis-Michael Lux <denis.lux@icloud.com>
 * @date    12.01.15
 */
namespace wigner
{
    // Wigner-d function
    auto wigner_d(const unsigned int& J, const int& M, const int& Mp, const double& beta) -> const double;
    
    // L^2 normalized Wiger-d function
    auto wigner_d_l2normalized(const unsigned int& J, const int& M, const int& Mp, const double& beta) -> const double;
    
    // Wigner-D function
    auto wigner_D(const unsigned int& J, const int& M, const int& Mp, const double& alpha, const double& beta, const double& gamma) -> const complex< double >;
    
    // L^2 normalized Wigner-D function
    auto wigner_D_l2normalized(const unsigned int& J, const int& M, const int& Mp, const double& alpha, const double& beta, const double& gamma) -> const complex< double >;
}

#endif