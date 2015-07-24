//
//  fn_rand_coeff.cpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 18.06.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_fn_rand_coeff_cpp
#define uzlmath_fn_rand_coeff_cpp

#include <uzlmath>

UZLMATH_BEGIN

/*!
 * @brief           filling Fourier coefficients container with random values
 * @details         Fills the fourier coefficients container with random values
 *                  in the given range.
 *
 * @param[in, out]  fc The fourier coefficients container
 * @param[in]       min The minimum value of the random number range
 * @param[in]       max The maximum value of the random number range
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            18.06.15
 *
 * @since           0.1.1
 */
auto rand(SOFTFourierCoefficients& fc, const double& min, const double& max) -> void
{
    if (min > max)
    {
        uzlmath_warning("min value is greater than max value in rand function for SOFTFourierCoefficients.");
        return;
    }
    
    // create timeval object
    struct timeval tv;
    
    // get current time in microseconds
    gettimeofday(&tv, NULL);
    
    // create seed
    unsigned long seed = 1000000L * tv.tv_sec + tv.tv_usec;
    
    // C++11 random numbers uniformly distributed
    std::default_random_engine e(seed);
    std::uniform_real_distribution< double > d(min, max);
    
    // iterate over degree
    for (int l = 0; l < fc.bandwidth; ++l)
    {
        // iterate over M order
        for (int M = -l; M <= l; ++M)
        {
            // iterate over M' order
            for (int Mp = -l; Mp <= l; ++Mp)
            {
                fc(l,M,Mp).re = d(e);
                fc(l,M,Mp).im = d(e);
            }
        }
    }
}

UZLMATH_END

#endif