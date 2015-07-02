//
//  forward_declarations.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 25.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_forward_declarations_hpp
#define uzlmath_forward_declarations_hpp

template< typename derived >            struct base;
template< typename T1, typename T2 >    class  glue;
template< typename eT >                 class  complex;

template< typename eT >                 class  matrix;
template< typename eT >                 class  matrix< complex< eT > >;

template< typename eT >                 class  vector;
template< typename eT >                 class  vector< complex< eT > >;

template< typename eT >                 class  grid3D;

template< typename eT >                 class  memory;

                                        class  factorial;
                                        class  stopwatch;

                                        struct SOFTFourierCoefficients;

/*!
 * @brief       Specifing the type a vector can have
 * @details     Indicator for the type of the current vector.
 *              A vector can either be a column vector or a
 *              row vector.
 *
 * @ingroup     vector
 */
enum class vec_type
{
    ROW,        //!< Represents a column vector of dimension \f$M\times 1\f$ where \f$M\in\mathbb{N}^+\f$
    COLUMN      //!< Represents a row vector of dimension \f$1\times M\f$ where \f$M\in\mathbb{N}^+\f$
};

/*!
 * @brief       Injection operator tokens for matrices
 * @details     This enum class contains tokens to manipulate 
 *              the injection operator for classes.
 */
enum class mat
{
    nextRow     //!< Represents a token indicating the next row for injection operators
};

#endif
