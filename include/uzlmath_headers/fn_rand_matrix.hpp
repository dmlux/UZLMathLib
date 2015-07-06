//
//  fn_rand_matrix.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 01.06.15.
//
//

#ifndef uzlmath_fn_rand_matrix_hpp
#define uzlmath_fn_rand_matrix_hpp

UZLMATH_BEGIN

/*!
 * @brief           Fills a given matrix with random integer values.
 * @details         The given vector gets filled with random integer values
 *                  in range of \f$[\min, \max]\f$.
 *
 * @param[in, out]  mat The matrix that is supposed to be filled with random
 *                  values.
 * @param[in]       max The max value for the random co-domain.
 * @param[in]       min The min value for the random co-domain.
 * @tparam          eT The element type of the matrix. The type can be any
 *                  number type.
 *
 * @since           0.1.1
 * @date            31.05.2015
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 *
 * @ingroup         matrix
 */
template< typename eT >
inline
auto randi(matrix< eT >& mat, const int& min, const int& max) -> typename uzl_void_real_num_only< eT >::result
{
    if (min > max)
    {
        printf("** uzlmath error: min value is greater than max value in rand function for SOFTFourierCoefficients. **");
        exit(EXIT_FAILURE);
    }
    
    // create timeval object
    struct timeval tv;
    
    // get current time in microseconds
    gettimeofday(&tv, NULL);
    
    // create seed
    unsigned long seed = 1000000L * tv.tv_sec + tv.tv_usec;
    
    // C++11 random numbers uniformly distributed
    std::default_random_engine e(seed);
    std::uniform_int_distribution< typename uzl_int_rand_dist_type< eT >::result > d(min, max);
    
    // fill vector with random integers
    size_t i, j;
    for (i = 0; i < mat.n_rows(); ++i)
    {
        for (j = 0; j < mat.n_cols(); ++j)
        {
            mat(i, j) = d(e);
        }
    }
}

/*!
 * @brief           Fills a given complex matrix with random integer values.
 * @details         The given vector gets filled with random integer values
 *                  in range of \f$[\min, \max]\f$.
 *
 * @param[in, out]  mat The complex matrix that is supposed to be filled with 
 *                  random values.
 * @param[in]       max The max value for the random co-domain.
 * @param[in]       min The min value for the random co-domain.
 * @tparam          eT The element type of the matrix. The type can be any
 *                  number type.
 *
 * @since           0.1.1
 * @date            31.05.2015
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 *
 * @ingroup         matrix
 */
template< typename eT >
inline
auto randi(matrix< complex< eT > >& mat, const int& min, const int& max) -> typename uzl_void_real_num_only< eT >::result
{
    if (min > max)
    {
        printf("** uzlmath error: min value is greater than max value in rand function for SOFTFourierCoefficients. **");
        exit(EXIT_FAILURE);
    }
    
    // create timeval object
    struct timeval tv;
    
    // get current time in microseconds
    gettimeofday(&tv, NULL);
    
    // create seed
    unsigned long seed = 1000000L * tv.tv_sec + tv.tv_usec;

    // C++11 random numbers uniformly distributed
    std::default_random_engine e(seed);
    std::uniform_int_distribution< typename uzl_int_rand_dist_type< eT >::result > d(min, max);
    
    // fill vector with random complex integers
    size_t i, j;
    for (i = 0; i < mat.n_rows(); ++i)
    {
        for (j = 0; j < mat.n_cols(); ++j)
        {
            mat(i, j).re = d(e);
            mat(i, j).im = d(e);
        }
    }
}

/*!
 * @brief           Fills a given real matrix with real random values.
 * @details         The given matrix gets filled with random real values in
 *                  range of \f$[\min, \max]\f$
 *
 * @param[in, out]  mat The real matrix that is supposed to be filled with random
 *                  values.
 * @param[in]       min The min value for the random co-domain.
 * @param[in]       max The max value for the random co-domain.
 * @tparam          eT The element type of the matrix. The type has to be a floating
 *                  point type (float, double or long double).
 *
 * @note            Only available for matrices of type float, double or long double!
 *
 * @since           0.1.1
 * @date            31.05.2015
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 *
 * @ingroup         matrix
 */
template< typename eT >
inline
auto rand(matrix< eT >& mat, const double& min, const double& max) -> typename uzl_void_real_only< eT >::result
{
    if (min > max)
    {
        printf("** uzlmath error: min value is greater than max value in rand function for SOFTFourierCoefficients. **");
        exit(EXIT_FAILURE);
    }
    
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
    size_t i, j;
    for (i = 0; i < mat.n_rows(); ++i)
    {
        for (j = 0; j < mat.n_cols(); ++j)
        {
            mat(i, j) = d(e);
        }
    }
}

/*!
 * @brief           Fills a given complex matrix with real random values.
 * @details         The given complex matrix gets filled with random real values in
 *                  range of \f$[\min, \max]\f$
 *
 * @param[in, out]  mat The complex matrix that is supposed to be filled with random
 *                  values.
 * @param[in]       min The min value for the random co-domain.
 * @param[in]       max The max value for the random co-domain.
 * @tparam          eT The element type of the matrix. The type has to be a floating
 *                  point type (float, double or long double).
 *
 * @note            Only available for matrices of type float, double or long double!
 *
 * @since           0.1.1
 * @date            31.05.2015
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 *
 * @ingroup         matrix
 */
template< typename eT >
inline
auto rand(matrix< complex< eT > >& mat, const double& min, const double& max) -> typename uzl_void_real_only< eT >::result
{
    if (min > max)
    {
        printf("** uzlmath error: min value is greater than max value in rand function for SOFTFourierCoefficients. **");
        exit(EXIT_FAILURE);
    }
    
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
    size_t i, j;
    for (i = 0; i < mat.n_rows(); ++i)
    {
        for (j = 0; j < mat.n_cols(); ++j)
        {
            mat(i, j).re = d(e);
            mat(i, j).im = d(e);
        }
    }
}

UZLMATH_END

#endif
