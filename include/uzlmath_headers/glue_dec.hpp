//
//  glue_dec.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 20.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_glue_dec_hpp
#define uzlmath_glue_dec_hpp

UZLMATH_BEGIN

/*!
 * @brief           Collection of struct and classes for the glue object
 * @defgroup        glue Glue class
 * @{
 */

/*!
 * @brief           The glue type represents two glued data structes for a special operation
 * @details         Combining two base objects and stores the references to them. Allows 
 *                  recursive definitions based on the fact that a glue object inherits from
 *                  the same base class than the data structures it should represent.
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @since           0.0.1
 */
template< typename T1, typename T2 >
class glue : public base< glue< T1, T2 > > // use static polymorphism
{
public:
    const T1& A;    //!< Constant reference to the first glue object
    const T2& B;    //!< Constant reference to the second glue object
    
    inline glue(const T1& in_A, const T2& in_B);
};

/*!
 * @brief           Helper struct to get depth of glue expression tree.
 * @details         The base case for the recursion struct.
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @since           0.1.1
 */
template< typename T1 >
struct depth_lhs
{
    static const int num = 0;   //!< Variable representing the depth of the BET
};

/*!
 * @brief           Helper struct to get depth of glue expression tree.
 * @details         The reurcive definition for the depth calculation. Since
 *                  expressions are getting evaluated from left to right the
 *                  The right object in the Glue type is always a plain object
 *                  and the left object can either be an other glue object of
 *                  next depth or in case of a base case a plain object.
 * 
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @since           0.1.1
 */
template< typename T1, typename T2 >
struct depth_lhs< glue< T1, T2 > >
{
    static const int num = 1 + depth_lhs<T1>::num;  // Try expanding the left node which is propably a glue type
};

/*!
 * @brief           Helper struct to extract all matrices from the binary expression
 *                  tree that represents the mathematical expression.
 * @details         This struct represents the base case of the BET.
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @since           0.1.1
 */
template< typename T1, typename eT >
struct mat_ptrs
{
    static const int num = 0;   //!< Variable representing the \f$N\f$-th element in expression
    
    inline static void get_ptrs(const matrix< eT >** ptrs, const T1& X);
};

/*!
 * @brief           Helper struct to extract all matrices from the binary expression
 *                  tree that represents the mathematical expression.
 * @details         Iterates recursively through the binary expression tree constructed
 *                  by glue objects and stores the pointers to the matrices in the given
 *                  matrix array.
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @since           0.1.1
 */
template< typename T1, typename T2, typename eT >
struct mat_ptrs< glue< T1, T2 >, eT >
{
    static const int num = 1 + mat_ptrs< T1, eT >::num; //!< Variable representing the \f$N\f$-th element in expression.
    
    inline static void get_ptrs(const matrix< eT >**in_ptrs, const glue< T1, T2 >& X);
};

/*!
 * @}
 */

UZLMATH_END

#endif
