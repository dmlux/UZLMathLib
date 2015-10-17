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
class matrix : public base< matrix< T > > // use static polymorphism
{
    size_t r_inj;       //!< Row injection index
    size_t c_inj;       //!< Column injection index
    
public:
    // ivars
    const size_t rows;  //!< Number of matrix rows
    const size_t cols;  //!< Number of matrix columns
    
    const T* mem;       //!< Matrix memory
    
    // methods
    inline                                         ~matrix();
    inline                                          matrix();
    
    inline                                          matrix(const size_t& m, const size_t& n);
    inline                                          matrix(const size_t& mn);
    inline                                          matrix(const size_t& m, const size_t& n, const T& initial);
    
    template< typename T1, typename T2 >
    inline                                          matrix(const glue< T1, T2 >& X);
    inline                                          matrix(const matrix< T >& A);
    inline                                          matrix(matrix< T >&& A);
    
    inline                                          matrix(const vector< T >& v);
    
    inline       glue< matrix< T >, matrix< T > >   operator+(const matrix< T >& A);
    inline       matrix< T >                        operator-(const matrix< T >& A);
    inline       matrix< T >                        operator*(const matrix< T >& A);
    inline       matrix< T >                        operator%(const matrix< T >& A);
    inline       matrix< T >                        operator/(const matrix< T >& A);
    inline       matrix< complex< T > >             operator+(const matrix< complex< T > >& A);
    inline       matrix< complex< T > >             operator-(const matrix< complex< T > >& A);
    inline       matrix< complex< T > >             operator*(const matrix< complex< T > >& A);
    inline       matrix< complex< T > >             operator%(const matrix< complex< T > >& A);
    inline       matrix< complex< T > >             operator/(const matrix< complex< T > >& A);
    
    inline       vector< T >                        operator*(const vector< T >& v);
    inline       vector< complex< T > >             operator*(const vector< complex< T > >& v);
    
    inline const matrix< T >&                       operator=(const matrix< T >& A);
    inline const matrix< T >&                       operator=(matrix< T >&& A);
    
    template< typename T1, typename T2 >
    inline const matrix< T >&                       operator=(const glue< T1, T2 >& X);
    
    inline       bool                               operator==(const matrix< T >& A);
    inline       bool                               operator!=(const matrix< T >& A);
    
    inline const matrix< T >&                       operator+=(const matrix< T >& A);
    inline const matrix< T >&                       operator-=(const matrix< T >& A);
    inline const matrix< T >&                       operator*=(const matrix< T >& A);
    inline const matrix< T >&                       operator%=(const matrix< T >& A);
    inline const matrix< T >&                       operator/=(const matrix< T >& A);
    
    inline const matrix< T >&                       operator*=(const vector< T >& v);
    
    inline       matrix< T >                        operator+();
    inline       matrix< T >                        operator-();
    
    inline       matrix< T >                        operator+(const T& rhs);
    inline       matrix< T >                        operator-(const T& rhs);
    inline       matrix< T >                        operator*(const T& rhs);
    inline       matrix< T >                        operator/(const T& rhs);
    inline       matrix< T >                        operator^(const unsigned int& exp);
    
    inline       matrix< complex< T > >             operator+(const complex< T >& rhs);
    inline       matrix< complex< T > >             operator-(const complex< T >& rhs);
    inline       matrix< complex< T > >             operator*(const complex< T >& rhs);
    inline       matrix< complex< T > >             operator/(const complex< T >& rhs);
    
    inline       matrix< T >&                       operator+=(const T& rhs);
    inline       matrix< T >&                       operator-=(const T& rhs);
    inline       matrix< T >&                       operator*=(const T& rhs);
    inline       matrix< T >&                       operator/=(const T& rhs);
    inline       matrix< T >&                       operator^=(const unsigned int& exp);
    
    inline       matrix< T >&                       operator<<(const T& val);
    inline       matrix< T >&                       operator<<(const mat& token);
    
    inline       T&                                 operator()(const size_t& i, const size_t& j);
    inline const T&                                 operator()(const size_t& i, const size_t& j) const;
    
    inline       void                               zeros();
    inline       void                               ones();
    inline       void                               eye();
    inline       void                               transpose();
    inline       void                               fill(const T& value);
    inline       void                               diag(const T& val);
    inline       void                               diag(const vector< T >& vec);
    inline       void                               reshape(const size_t& new_rows, const size_t& new_cols);
    
    inline const double                             determinant();
};

template< typename S > std::ostream& operator<<(std::ostream& o, const matrix< S >& A);

/*!
 * @}
 */

UZLMATH_END

#endif
