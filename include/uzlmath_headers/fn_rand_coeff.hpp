//
//  fn_rand_coeff.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 06.06.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_fn_rand_coeff_hpp
#define uzlmath_fn_rand_coeff_hpp

UZLMATH_BEGIN

// filling SOFTFourierCoefficients container with randoms in range
auto rand(SOFTFourierCoefficients& fc, const double& min, const double& max) -> void;

UZLMATH_END

#endif
