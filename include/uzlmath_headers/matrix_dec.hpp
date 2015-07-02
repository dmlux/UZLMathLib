//
//  matrix_tpl_dec.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 13.01.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_matrix_dec_hpp
#define uzlmath_matrix_dec_hpp

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
 * @tparam      eT An element type which represents a number that provides all common
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
template< typename eT >
class matrix : public base< matrix< eT > > // use static polymorphism
{
    size_t rows;                //!< Number of matrix rows
    size_t cols;                //!< Number of matrix columns
    
    size_t r_inj;               //!< Row injection index
    size_t c_inj;               //!< Column injection index
    
    eT* mem;                    //!< Matrix memory
    
    // Friend classes to grant member access
    template< typename F > friend class matrix;
    template< typename F > friend class vector;
    
public:
    inline                                         ~matrix();
    inline                                          matrix();
    
    inline                                          matrix(const size_t& m, const size_t& n);
    inline                                          matrix(const size_t& mn);
    inline                                          matrix(const size_t& m, const size_t& n, const eT& initial);
    
    template< typename T1, typename T2 >
    inline                                          matrix(const glue< T1, T2 >& X);
    inline                                          matrix(const matrix< eT >& A);
    inline                                          matrix(matrix< eT >&& A);
    
    inline                                          matrix(const vector< eT >& v);
    
    inline       glue< matrix< eT >, matrix< eT > > operator+(const matrix< eT >& A);
    inline       matrix< eT >                       operator-(const matrix< eT >& A);
    inline       matrix< eT >                       operator*(const matrix< eT >& A);
    inline       matrix< eT >                       operator%(const matrix< eT >& A);
    inline       matrix< eT >                       operator/(const matrix< eT >& A);
    inline       matrix< complex< eT > >            operator+(const matrix< complex< eT > >& A);
    inline       matrix< complex< eT > >            operator-(const matrix< complex< eT > >& A);
    inline       matrix< complex< eT > >            operator*(const matrix< complex< eT > >& A);
    inline       matrix< complex< eT > >            operator%(const matrix< complex< eT > >& A);
    inline       matrix< complex< eT > >            operator/(const matrix< complex< eT > >& A);
    
    inline       vector< eT >                       operator*(const vector< eT >& v);
    inline       vector< complex< eT > >            operator*(const vector< complex< eT > >& v);
    
    inline const matrix< eT >&                      operator=(const matrix< eT >& A);
    inline const matrix< eT >&                      operator=(matrix< eT >&& A);
    
    template< typename T1, typename T2 >
    inline const matrix< eT >&                      operator=(const glue< T1, T2 >& X);
    
    inline       bool                               operator==(const matrix< eT >& A);
    inline       bool                               operator!=(const matrix< eT >& A);
    
    inline const matrix< eT >&                      operator+=(const matrix< eT >& A);
    inline const matrix< eT >&                      operator-=(const matrix< eT >& A);
    inline const matrix< eT >&                      operator*=(const matrix< eT >& A);
    inline const matrix< eT >&                      operator%=(const matrix< eT >& A);
    inline const matrix< eT >&                      operator/=(const matrix< eT >& A);
    
    inline const matrix< eT >&                      operator*=(const vector< eT >& v);
    
    inline       matrix< eT >                       operator+();
    inline       matrix< eT >                       operator-();
    
    inline       matrix< eT >                       operator+(const eT& rhs);
    inline       matrix< eT >                       operator-(const eT& rhs);
    inline       matrix< eT >                       operator*(const eT& rhs);
    inline       matrix< eT >                       operator/(const eT& rhs);
    inline       matrix< eT >                       operator^(const unsigned int& exp);
    
    inline       matrix< complex< eT > >            operator+(const complex< eT >& rhs);
    inline       matrix< complex< eT > >            operator-(const complex< eT >& rhs);
    inline       matrix< complex< eT > >            operator*(const complex< eT >& rhs);
    inline       matrix< complex< eT > >            operator/(const complex< eT >& rhs);
    
    inline       matrix< eT >&                      operator+=(const eT& rhs);
    inline       matrix< eT >&                      operator-=(const eT& rhs);
    inline       matrix< eT >&                      operator*=(const eT& rhs);
    inline       matrix< eT >&                      operator/=(const eT& rhs);
    inline       matrix< eT >&                      operator^=(const unsigned int& exp);
    
    inline       matrix< eT >&                      operator<<(const eT& val);
    inline       matrix< eT >&                      operator<<(const mat& T);
    
    inline       eT&                                operator()(const size_t& i, const size_t& j);
    inline constexpr eT&                            operator()(const size_t& i, const size_t& j) const;
    
    inline       void                               zeros();
    inline       void                               ones();
    inline       void                               eye();
    inline       void                               transpose();
    inline       void                               fill(const eT& value);
    inline       void                               diag(const eT& val);
    inline       void                               diag(const vector< eT >& vec);
    
    inline constexpr size_t                         n_rows() const;
    inline constexpr size_t                         n_cols() const;
    inline       eT*                                memptr();
    inline const eT*                                memptr() const;
    
    inline const double                             determinant();
};

template< typename S > std::ostream& operator<<(std::ostream& o, const matrix< S >& A);

/*!
 * @}
 */

#endif
