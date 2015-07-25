//
//  fn_obj2file_soft_coef.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 16.06.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_fn_obj2file_soft_coef_hpp
#define UZLMathLib_fn_obj2file_soft_coef_hpp

UZLMATH_BEGIN

// Saves a SOFTFourierCoefficients container on disk
auto obj2file(const SOFTFourierCoefficients& coef, const std::string& fileName) -> void;

UZLMATH_END

#endif
