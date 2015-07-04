//
//  fn_rand_vector.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 31.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_fn_rand_vector_hpp
#define uzlmath_fn_rand_vector_hpp

/*!
 * @brief           Fills a given vector with random integer values.
 * @details         The given vector gets filled with random integer values
 *                  in range of \f$[\min, \max]\f$.
 * 
 * @param[in, out]  vec The vector that is supposed to be filled with random
 *                  values.
 * @param[in]       max The max value for the random co-domain.
 * @param[in]       min The min value for the random co-domain.
 *
 * @ingroup         vector
 */
template< typename eT >
inline
auto randi(vector< eT >& vec, const int& min, const int& max) -> typename uzl_void_real_num_only< eT >::result
{
    // create timeval object
    struct timeval tv;
    
    // get current time in microseconds
    gettimeofday(&tv, NULL);
    
    // create seed
    unsigned long seed = 1000000L * tv.tv_sec + tv.tv_usec;
    
    // C++11 random numbers uniformly distributed
    std::default_random_engine e(seed);
    std::uniform_int_distribution< typename uzl_int_rand_dist_type< eT >::result > d(min, max);
    
    // fill vector with randoms
    size_t i;
    for (i = 0; i < vec.n_elements(); ++i)
    {
        vec[i] = d(e);
    }
}

template< typename eT >
inline
auto randi(vector< complex< eT > >& vec, const int& min, const int& max) -> typename uzl_void_real_num_only< eT >::result
{
    // create timeval object
    struct timeval tv;
    
    // get current time in microseconds
    gettimeofday(&tv, NULL);
    
    // create seed
    unsigned long seed = 1000000L * tv.tv_sec + tv.tv_usec;
    
    // C++11 random numbers uniformly distributed
    std::default_random_engine e(seed);
    std::uniform_int_distribution< typename uzl_int_rand_dist_type< eT >::result > d(min, max);
    
    // fill vector with randoms
    size_t i;
    for (i = 0; i < vec.n_elements(); ++i)
    {
        vec[i].re = d(e);
        vec[i].im = d(e);
    }
}


/*!
 * @brief           Fills a given real vector with real random values.
 * @details         The given vector gets filled with random real values in
 *                  range of \f$[\min, \max]\f$
 *
 * @param[in, out]  vec The real vector that is supposed to be filled with random
 *                  values.
 * @param[in]       min The min value for the random co-domain.
 * @param[in]       max The max value for the random co-domain.
 * @tparam          eT The element type of the vector. The type has to be a floating 
 *                  point type (float, double or long double).
 *
 * @note            Only available for vectors of type float, double or long double!
 *
 * @ingroup         vector
 */
template< typename eT >
auto rand(vector< eT >& vec, const double& min, const double& max) -> typename uzl_void_real_only< eT >::result
{
    // create timeval object
    struct timeval tv;
    
    // get current time in microseconds
    gettimeofday(&tv, NULL);
    
    // create seed
    unsigned long seed = 1000000L * tv.tv_sec + tv.tv_usec;
    
    // C++11 random numbers uniformly distributed
    std::default_random_engine e(seed);
    std::uniform_real_distribution< eT > d(min, max);
    
    // fill vector with randoms
    size_t i;
    for (i = 0; i < vec.n_elements(); ++i)
    {
        vec[i] = d(e);
    }
}

/*!
 * @brief           Fills a given complex vector with real random values.
 * @details         The given vector gets filled with random complex values in
 *                  range of \f$[\min, \max]\f$
 *
 * @param[in, out]  vec The complex vector that is supposed to be filled with random
 *                  values.
 * @param[in]       min The min value for the random co-domain.
 * @param[in]       max The max value for the random co-domain.
 * @tparam          eT The element type of the vector. The type has to be a floating
 *                  point type (float, double or long double).
 *
 * @note            Only available for vectors of type float, double or long double!
 *
 * @ingroup         vector
 */
template< typename eT >
auto rand(vector< complex< eT > >& vec, const double& min, const double& max) -> typename uzl_void_real_only< eT >::result
{
    // create timeval object
    struct timeval tv;
    
    // get current time in microseconds
    gettimeofday(&tv, NULL);
    
    // create seed
    unsigned long seed = 1000000L * tv.tv_sec + tv.tv_usec;
    
    // C++11 random numbers uniformly distributed
    std::default_random_engine e(seed);
    std::uniform_real_distribution< eT > d(min, max);
    
    // fill vector with randoms
    size_t i;
    for (i = 0; i < vec.n_elements(); ++i)
    {
        vec[i].re = d(e);
        vec[i].im = d(e);
    }
}

#endif
