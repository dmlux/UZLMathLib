//
//  fn_fliplr.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 16.06.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_fn_fliplr_hpp
#define uzlmath_fn_fliplr_hpp

/*!
 * @brief           Flips a given matrix by mirroring the columns.
 * @details         Flips a matrix so that 
 *                  \f$M = (v_0,v_1,v_2,\dots,v_{n-2}, v_{n-1}, v_{n})\f$
 *                  becomes
 *                  \f$M = (v_n,v_{n-1},v_{n-2},\dots,v_2,v_1,v_0)\f$
 *
 * @param[in,out]   mat The matrix that is supposed to be flipped
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            16.06.15
 *
 * @since           0.1.1
 *
 * @ingroup         matrix
 */
template< typename eT >
inline
auto fliplr(matrix< eT >& mat) -> typename uzl_void_real_num_only< eT >::result
{
    // define indices
    size_t j, k;
    
    // iterate over half of columns
    for (j = 0; j < mat.n_cols() / 2; ++j)
    {
        // iterate over rows: Since the matrix is in column major order
        // in memory this is important to get maximum speed!
        for (k = 0; k < mat.n_rows(); ++k)
        {
            // temporary element copy
            eT tmp                       = mat(k, j);
            // swap elements
            mat(k, j)                    = mat(k, mat.n_cols() - j - 1);
            mat(k, mat.n_cols() - j - 1) = tmp;
        }
    }
}

/*!
 * @brief           Flips a given matrix by inverting the column order and negating
 *                  each element in every second even row.
 * @details         Flips a matrix so that
 *                  \f$M = (v_0,v_1,v_2,\dots,v_{n-2}, v_{n-1}, v_{n})\f$
 *                  becomes
 *                  \f$M = (v_n,v_{n-1},v_{n-2},\dots,v_2,v_1,v_0)\f$
 *                  with negated elements in every even row.
 *
 * @param[in,out]   mat The matrix that is supposed to be flipped
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            16.06.15
 *
 * @since           0.1.1
 *
 * @ingroup         matrix
 */
template< typename eT >
inline
auto fliplr_ne2nderow(matrix< eT >& mat) -> typename uzl_void_real_num_only<eT>::result
{
    // define indices
    size_t j, k;
    
    // iterate over half of columns
    for (j = 0; j < mat.n_cols() / 2; ++j)
    {
        // iterate over rows: Since the matrix is in column major order
        // in memory this is important to get maximum speed!
        for (k = 0; k < mat.n_rows(); ++k)
        {
            // temporary element copy
            eT tmp                       = mat(k, j);
            
            // swap elements
            if (k & 1)
            {
                mat(k, j)                    = mat(k, mat.n_cols() - j - 1);
                mat(k, mat.n_cols() - j - 1) = tmp;
            }
            else
            {
                mat(k, j)                    =  -mat(k, mat.n_cols() - j - 1);
                mat(k, mat.n_cols() - j - 1) = -tmp;
            }
        }
    }
}

/*!
 * @brief           Flips a given matrix by inverting the column order and negating
 *                  each element in every second odd row.
 * @details         Flips a matrix so that
 *                  \f$M = (v_0,v_1,v_2,\dots,v_{n-2}, v_{n-1}, v_{n})\f$
 *                  becomes
 *                  \f$M = (v_n,v_{n-1},v_{n-2},\dots,v_2,v_1,v_0)\f$
 *                  with negated element in every odd row.
 *
 * @param[in,out]   mat The matrix that is supposed to be flipped
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            16.06.15
 *
 * @since           0.1.1
 *
 * @ingroup         matrix
 */
template< typename eT >
inline
auto fliplr_ne2ndorow(matrix< eT >& mat) -> typename uzl_void_real_num_only<eT>::result
{
    // define indices
    size_t j, k;
    
    // iterate over half of columns
    for (j = 0; j < mat.n_cols() / 2; ++j)
    {
        // iterate over rows: Since the matrix is in column major order
        // in memory this is important to get maximum speed!
        for (k = 0; k < mat.n_rows(); ++k)
        {
            // temporary element copy
            eT tmp                       = mat(k, j);
            
            // swap elements
            if (k & 1)
            {
                mat(k, j)                    =  -mat(k, mat.n_cols() - j - 1);
                mat(k, mat.n_cols() - j - 1) = -tmp;
            }
            else
            {
                mat(k, j)                    = mat(k, mat.n_cols() - j - 1);
                mat(k, mat.n_cols() - j - 1) = tmp;
            }
        }
    }
}

/*!
 * @brief           Flips a given complex matrix by mirroring the columns.
 * @details         Flips a complex matrix so that
 *                  \f$M = (v_0,v_1,v_2,\dots,v_{n-2}, v_{n-1}, v_{n})\f$
 *                  becomes
 *                  \f$M = (v_n,v_{n-1},v_{n-2},\dots,v_2,v_1,v_0)\f$
 *
 * @param[in,out]   mat The complex matrix that is supposed to be flipped
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            16.06.15
 *
 * @since           0.1.1
 *
 * @ingroup         matrix
 */
template< typename eT >
inline
auto fliplr(matrix< complex< eT > >& mat) -> typename uzl_void_real_num_only< eT >::result
{
    // define indices
    size_t j, k;
    
    // iterate over half of columns
    for (j = 0; j < mat.n_cols() / 2; ++j)
    {
        // iterate over rows: Since the matrix is in column major order
        // in memory this is important to get maximum speed!
        for (k = 0; k < mat.n_rows(); ++k)
        {
            // temporary element copy
            complex< eT > tmp            = mat(k, j);
            // swap elements
            mat(k, j)                    = mat(k, mat.n_cols() - j - 1);
            mat(k, mat.n_cols() - j - 1) = tmp;
        }
    }
}

#endif
