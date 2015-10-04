//
//  glue_def.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 23.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_glue_def_hpp
#define UZLMathLib_glue_def_hpp

UZLMATH_BEGIN

/*!
 * @brief           The glue type default constructor.
 * @details         Builds the glue object by just applying both elements of the expression to
 *                  the glue element ivars.
 */
template< typename T1, typename T2 >
inline
glue< T1, T2 >::glue(const T1& in_A, const T2& in_B)
    : A(in_A)
    , B(in_B)
{}



/*!
 * @brief           Recursive function that stores the matrix from the current binary expression.
 * @details         This function represents the base case where only two matrices are existing in
 *                  In the current binary expression.
 */
template< typename T1, typename T >
inline
void mat_ptrs< T1, T >::get_ptrs(const matrix< T >** ptrs, const T1& X)
{
    ptrs[0] = reinterpret_cast< const matrix< T >* >(&X);
}



/*!
 * @brief           Recursive function that stores the matrix from the current binary expression.
 * @details         This function traverse the BET and stores on each level the matrix element and
 *                  calls itself with the other binary expression.
 * 
 * @param[in]       in_ptrs The matrix pointer array to store the matrix reference.
 * @param[in]       X The current expression in the tree
 */
template< typename T1, typename T2, typename T >
inline
void mat_ptrs< glue< T1, T2 >, T >::get_ptrs(const matrix< T >**in_ptrs, const glue< T1, T2 >& X)
{
    // traverse left nodes
    mat_ptrs< T1, T >::get_ptrs(in_ptrs, X.A);
    
    // get address of the matrix on the right
    in_ptrs[num] = reinterpret_cast< const matrix< T >* >(&X.B);
}

UZLMATH_END

#endif
