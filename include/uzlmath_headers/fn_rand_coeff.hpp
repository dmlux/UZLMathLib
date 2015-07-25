//
//  fn_rand_coeff.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 06.06.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_fn_rand_coeff_hpp
#define UZLMathLib_fn_rand_coeff_hpp

UZLMATH_BEGIN

// filling SOFTFourierCoefficients container with randoms in range
auto rand(SOFTFourierCoefficients& fc, const double& min, const double& max) -> void;

UZLMATH_END

#endif
