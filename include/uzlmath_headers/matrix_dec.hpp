//
//  matrix_dec.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 13.01.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_matrix_dec_hpp
#define UZLMathLib_matrix_dec_hpp

UZLMATH_BEGIN

/*!
 * @brief       Collection of classes and functions for matrices for mathematical
 *              purposes.
 * @defgroup    matrix Matrix
 * @{
 */

/*!
 * @brief       Matrix class with common operations.
 * @details     The matrix class uses internally an dynamic allocated array
 *              to store values. The memory is organized in  column-major
 *              order to provide interchangability with the BLAS and LAPACK
 *              Fortran libraries.
 *
 * @tparam      T An element type which represents a number that provides all common
 *              mathmatical operations.
 *
 * @since       0.0.1
 *
 * @todo        Implement matrix transposition with flags to change setting
 *              and getting behavior of the ()-operator. This will save some
 *              extra computation
 * @todo        Implement injection functionality and token for the <<-Operator.
 *
 * @author      Denis-Michael Lux <denis.lux@icloud.com>
 * @date        12.01.15
 */
template< typename T >
class
matrix< T, typename if_true< is_num_type< T >::value >::type > : public base< matrix< T > > // use static polymorphism
{
    size_t r_inj;           //!< Row injection index
    size_t c_inj;           //!< Column injection index
    
public:
    // type declarations
    typedef T pod_type;     //!< The POD type of matrix elements
    
    // ivars
    const size_t rows;      //!< Number of matrix rows
    const size_t cols;      //!< Number of matrix columns
    
    const pod_type* mem;    //!< Matrix memory
    
    // methods
    inline                                                     ~matrix();
    inline                                                      matrix();
    
    inline                                                      matrix(const size_t& m, const size_t& n);
    inline                                                      matrix(const size_t& mn);
    inline                                                      matrix(const size_t& m, const size_t& n, const pod_type& initial);
    
    template< typename T1, typename T2 >
    inline                                                      matrix(const glue< T1, T2 >& X);
    inline                                                      matrix(const matrix< pod_type >& A);
    inline                                                      matrix(matrix< pod_type >&& A);
    
    inline                                                      matrix(const vector< pod_type >& v);
    
    inline       glue< matrix< pod_type >, matrix< pod_type > > operator+(const matrix< pod_type >& A);
    inline       matrix< pod_type >                             operator-(const matrix< pod_type >& A);
    inline       matrix< pod_type >                             operator*(const matrix< pod_type >& A);
    inline       matrix< pod_type >                             operator%(const matrix< pod_type >& A);
    inline       matrix< pod_type >                             operator/(const matrix< pod_type >& A);
    inline       matrix< complex< pod_type > >                  operator+(const matrix< complex< pod_type > >& A);
    inline       matrix< complex< pod_type > >                  operator-(const matrix< complex< pod_type > >& A);
    inline       matrix< complex< pod_type > >                  operator*(const matrix< complex< pod_type > >& A);
    inline       matrix< complex< pod_type > >                  operator%(const matrix< complex< pod_type > >& A);
    inline       matrix< complex< pod_type > >                  operator/(const matrix< complex< pod_type > >& A);
    
    inline       vector< pod_type >                             operator*(const vector< pod_type >& v);
    inline       vector< complex< pod_type > >                  operator*(const vector< complex< pod_type > >& v);
    
    inline const matrix< pod_type >&                            operator=(const matrix< pod_type >& A);
    inline const matrix< pod_type >&                            operator=(matrix< pod_type >&& A);
    
    template< typename T1, typename T2 >
    inline const matrix< pod_type >&                            operator=(const glue< T1, T2 >& X);
    
    inline       bool                                           operator==(const matrix< pod_type >& A);
    inline       bool                                           operator!=(const matrix< pod_type >& A);
    
    inline const matrix< pod_type >&                            operator+=(const matrix< pod_type >& A);
    inline const matrix< pod_type >&                            operator-=(const matrix< pod_type >& A);
    inline const matrix< pod_type >&                            operator*=(const matrix< pod_type >& A);
    inline const matrix< pod_type >&                            operator%=(const matrix< pod_type >& A);
    inline const matrix< pod_type >&                            operator/=(const matrix< pod_type >& A);
    
    inline const matrix< pod_type >&                            operator*=(const vector< pod_type >& v);
    
    inline       matrix< pod_type >                             operator+();
    inline       matrix< pod_type >                             operator-();
    
    inline       matrix< pod_type >                             operator+(const pod_type& rhs);
    inline       matrix< pod_type >                             operator-(const pod_type& rhs);
    inline       matrix< pod_type >                             operator*(const pod_type& rhs);
    inline       matrix< pod_type >                             operator/(const pod_type& rhs);
    inline       matrix< pod_type >                             operator^(const unsigned int& exp);
    
    inline       matrix< complex< pod_type > >                  operator+(const complex< pod_type >& rhs);
    inline       matrix< complex< pod_type > >                  operator-(const complex< pod_type >& rhs);
    inline       matrix< complex< pod_type > >                  operator*(const complex< pod_type >& rhs);
    inline       matrix< complex< pod_type > >                  operator/(const complex< pod_type >& rhs);
    
    inline       matrix< pod_type >&                            operator+=(const pod_type& rhs);
    inline       matrix< pod_type >&                            operator-=(const pod_type& rhs);
    inline       matrix< pod_type >&                            operator*=(const pod_type& rhs);
    inline       matrix< pod_type >&                            operator/=(const pod_type& rhs);
    inline       matrix< pod_type >&                            operator^=(const unsigned int& exp);
    
    inline       matrix< pod_type >&                            operator<<(const pod_type& val);
    inline       matrix< pod_type >&                            operator<<(const mat& token);
    
    inline       pod_type&                                      operator()(const size_t& i, const size_t& j);
    inline const pod_type&                                      operator()(const size_t& i, const size_t& j) const;
    
    inline       void                                           zeros();
    inline       void                                           ones();
    inline       void                                           eye();
    inline       void                                           transpose();
    inline       void                                           fill(const pod_type& value);
    inline       void                                           diag(const pod_type& val);
    inline       void                                           diag(const vector< pod_type >& vec);
    inline       void                                           reshape(const size_t& new_rows, const size_t& new_cols);
    
    inline       double                                         determinant();
};

template< typename S > std::ostream& operator<<(std::ostream& o, const matrix< S >& A);

/*!
 * @}
 */

UZLMATH_END

#endif
