//
//  fn_spharmonics.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 15.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_fn_spharmonics_hpp
#define uzlmath_fn_spharmonics_hpp

/*!
 * @brief       **Spherical harmonics** are a series of special functions defined on the 
 *              surface of a sphere used to solve some kinds of differential equations.
 *
 * @since       0.0.1
 *
 * @author      Denis-Michael Lux <denis.lux@icloud.com>
 * @date        21.04.15
 */
namespace spharmonics
{
    // Forward fast Fourier transform on SO(3)
    auto SOFT(grid3D< complex< double > > sample, SOFTFourierCoefficients& fc) -> void;
    
    // Inverse fast Fourier transform on SO(3)
    auto inverseSOFT(const SOFTFourierCoefficients& fc, grid3D< complex< double > >& synthesis) -> void;
}

#endif
