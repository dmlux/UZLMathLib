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

//template< typename pod_type, typename T >
//void test_randctx(const randctx< pod_type, T > ctx)
//{
//    if (is_normal_distribution< T >::value == true)
//    {
//        printf("\nOK war ne Normalverteilung...\n");
//    }
//    else if (is_uniform_distribution< T >::value == true)
//    {
//        printf("\nOK war ne diskrete Gleichverteilung...\n");
//    }
//    else
//    {
//        printf("\nHmmm... war keine bekannte Verteilung...\n");
//    }
//}
//
//auto rand(SOFTFourierCoefficients& fc, const double& min, const double& max) -> void
//{
//    if (min > max)
//    {
//        uzlmath_warning("%s", "min value is greater than max value in rand function for SOFTFourierCoefficients.");
//        return;
//    }
//    
//    // create timeval object
//    struct timeval tv;
//    
//    // get current time in microseconds
//    gettimeofday(&tv, NULL);
//    
//    // create seed
//    unsigned long seed = 1000000L * tv.tv_sec + tv.tv_usec;
//    
//    // C++11 random numbers uniformly distributed with marsenne twister
//    std::mt19937_64 e(seed);
//    std::uniform_real_distribution< double > d(min, max);
//    
//    // iterate over degree
//    for (int l = 0; l < fc.bandwidth; ++l)
//    {
//        // iterate over M order
//        for (int M = -l; M <= l; ++M)
//        {
//            // iterate over M' order
//            for (int Mp = -l; Mp <= l; ++Mp)
//            {
//                fc(l,M,Mp).re = d(e);
//                fc(l,M,Mp).im = d(e);
//            }
//        }
//    }
//}

// filling SOFTFourierCoefficients container with randoms in range
auto rand(SOFTFourierCoefficients& fc, const double& min, const double& max) -> void;

template< typename pod_type, typename distribution >
inline
auto rand(SOFTFourierCoefficients& fc, const randctx< pod_type, distribution >& ctx) -> typename uzl_double_only< pod_type >::result
{
    
}

UZLMATH_END

#endif
