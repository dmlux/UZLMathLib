//
//  fn_flipud.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 03.07.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_fn_flipud_hpp
#define UZLMathLib_fn_flipud_hpp

UZLMATH_BEGIN

/*!
 * @brief           Flips a given matrix inverting the row order.
 * @details         Flips a matrix so that
 *                  \f$M = \begin{pmatrix}v_0\\ v_1\\ v_2\\ \dots\\ v_{n-2}\\ v_{n-1}\\ v_{n}\end{pmatrix}\f$
 *                  becomes
 *                  \f$M = \begin{pmatrix}v_n\\ v_{n-1}\\ v_{n-2}\\ \dots\\ v_2\\ v_1\\ v_0\end{pmatrix}\f$
 *
 * @param[in,out]   mat The matrix that is supposed to be flipped
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            03.07.15
 *
 * @since           0.1.1
 *
 * @ingroup         matrix
 */
template< typename eT >
inline
auto flipud(matrix< eT >& mat) -> typename uzl_void_real_num_only< eT >::result
{
    // define indices
    size_t j, k;
    
    // iterate over cols
    for (k = 0; k < mat.cols; ++k)
    {
        // iterate over half of rows: Since the matrix is in column major order
        // in memory this is important to get maximum speed!
        for (j = 0; j < mat.rows / 2; ++j)
        {
            // temporary element copy
            eT tmp = mat(j, k);
            // swap elements
            mat(j, k)                = mat(mat.rows - j - 1, k);
            mat(mat.rows - j - 1, k) = tmp;
        }
    }
}

/*!
 * @brief           Flips a given matrix by inverting the row order and negating
 *                  each element in every second even column.
 * @details         Flips a matrix so that
 *                  \f$M = (v_0,v_1,v_2,\dots,v_{n-2}, v_{n-1}, v_{n})\f$
 *                  becomes
 *                  \f$M = (v_n,v_{n-1},v_{n-2},\dots,v_2,v_1,v_0)\f$
 *                  with negated elements in every second even column.
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
auto flipud_ne2ndecol(matrix< eT >& mat) -> typename uzl_void_real_num_only< eT >::result
{
    // define indices
    size_t j, k;
    
    // iterate over cols
    for (k = 0; k < mat.cols; ++k)
    {
        // iterate over half of rows: Since the matrix is in column major order
        // in memory this is important to get maximum speed!
        for (j = 0; j < mat.rows / 2; ++j)
        {
            // temporary element copy
            eT tmp = mat(j, k);
            
            // swap elements
            if (k & 1)
            {
                mat(j, k)                = mat(mat.rows - j - 1, k);
                mat(mat.rows - j - 1, k) = tmp;
            }
            else
            {
                mat(j, k)                = -mat(mat.rows - j - 1, k);
                mat(mat.rows - j - 1, k) = -tmp;
            }
        }
    }
}

/*!
 * @brief           Flips a given matrix by inverting the row order and negating
 *                  each element in every second odd column.
 * @details         Flips a matrix so that
 *                  \f$M = (v_0,v_1,v_2,\dots,v_{n-2}, v_{n-1}, v_{n})\f$
 *                  becomes
 *                  \f$M = (v_n,v_{n-1},v_{n-2},\dots,v_2,v_1,v_0)\f$
 *                  with negated elements in every second odd column.
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
auto flipud_ne2ndocol(matrix< eT >& mat) -> typename uzl_void_real_num_only< eT >::result
{
    // define indices
    size_t j, k;
    
    // iterate over cols
    for (k = 0; k < mat.cols; ++k)
    {
        // iterate over half of rows: Since the matrix is in column major order
        // in memory this is important to get maximum speed!
        for (j = 0; j < mat.rows / 2; ++j)
        {
            // temporary element copy
            eT tmp = mat(j, k);
            
            // swap elements
            if (k & 1)
            {
                mat(j, k)                = -mat(mat.rows - j - 1, k);
                mat(mat.rows - j - 1, k) = -tmp;
            }
            else
            {
                mat(j, k)                = mat(mat.rows - j - 1, k);
                mat(mat.rows - j - 1, k) = tmp;
            }
        }
    }
}

/*!
 * @brief           Flips a given complex matrix inverting the row order.
 * @details         Flips a complex matrix so that
 *                  \f$M = \begin{pmatrix}v_0\\ v_1\\ v_2\\ \dots\\ v_{n-2}\\ v_{n-1}\\ v_{n}\end{pmatrix}\f$
 *                  becomes
 *                  \f$M = \begin{pmatrix}v_n\\ v_{n-1}\\ v_{n-2}\\ \dots\\ v_2\\ v_1\\ v_0\end{pmatrix}\f$
 *
 * @param[in,out]   mat The complex matrix that is supposed to be flipped
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            03.07.15
 *
 * @since           0.1.1
 *
 * @ingroup         matrix
 */
template< typename eT >
inline
auto flipud(matrix< complex< eT > >& mat) -> typename uzl_void_real_num_only< eT >::result
{
    // define indices
    size_t j, k;
    
    // iterate over half of columns
    for (k = 0; k < mat.n_cols(); ++k)
    {
        // iterate over half of rows: Since the matrix is in column major order
        // in memory this is important to get maximum speed!
        for (j = 0; j < mat.n_rows() / 2; ++j)
        {
            // temporary element copy
            complex< eT > tmp = mat(j, k);
            // swap elements
            mat(j, k)                    = mat(mat.n_rows() - j - 1, k);
            mat(mat.n_rows() - j - 1, k) = tmp;
        }
    }
}

UZLMATH_END

#endif
