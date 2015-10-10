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

// filling DSOFTFourierCoefficients container with randoms in range
void rand(DSOFTFourierCoefficients& fc, const double& min, const double& max);

template< typename pod_type, typename distribution >
inline
typename uzl_uniform_real_dist_only< distribution >::result rand(DSOFTFourierCoefficients& fc, const randctx< pod_type, distribution >& ctx)
{
    // cast to correct underlying type
    uniform_real_distribution< pod_type > dist = ctx.get_ref();
    
    // get min and max values
    pod_type min = dist.min < dist.max ? dist.min : dist.max;
    pod_type max = dist.min < dist.max ? dist.max : dist.min;
    
    // create timeval objects
    struct timeval tv;
    
    // get current time in microseconds
    gettimeofday(&tv, NULL);
    
    // create seed
    unsigned long seed = 1000000L * tv.tv_sec + tv.tv_usec;
    
    // C++11 random numbers uniformly distributed
    if ( dist.engine == random_engine::DEFAULT )
    {
        std::mt19937_64 e(seed);
        std::uniform_real_distribution< pod_type > d(min, max);
        
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
    
    else if (dist.engine == random_engine::MERSENNE_TWISTER)
    {
        std::mt19937 e(seed);
        std::uniform_real_distribution< pod_type > d(min, max);
        
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
}

UZLMATH_END

#endif
