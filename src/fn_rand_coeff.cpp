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

namespace uzlmath
{
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
    auto rand_coef(SOFTFourierCoefficients& fc, const double& min, const double& max) -> void
    {
        // create timeval object
        struct timeval tv;
        
        // get current time in microseconds
        gettimeofday(&tv, NULL);
        
        // create seed
        unsigned long seed = 1000000 * tv.tv_sec + tv.tv_usec;
        
        // seed the random number generator
        srand(seed);
        
        // drop first seed
        rand();
        
        // iterate over degree
        for (int l = 0; l < fc.bandwidth(); ++l)
        {
            // iterate over M order
            for (int M = -l; M <= l; ++M)
            {
                // iterate over M' order
                for (int Mp = -l; Mp <= l; ++Mp)
                {
                    fc(l,M,Mp).re = (static_cast< double >(rand()) / RAND_MAX) * (max - min) + min;
                    fc(l,M,Mp).im = (static_cast< double >(rand()) / RAND_MAX) * (max - min) + min;
                }
            }
        }
    }
}

#endif