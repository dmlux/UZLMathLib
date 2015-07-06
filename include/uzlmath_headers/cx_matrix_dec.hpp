//
//  cx_matrix_dec.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 13.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_cx_matrix_dec_hpp
#define uzlmath_cx_matrix_dec_hpp

UZLMATH_BEGIN

/*!
 * @brief       Matrix specialization for complex data.
 * @details     A matrix that contains only complex values on each
 *              entry.  The matrix class uses internally an dynamic allocated array
 *              to store values. The memory is organized in  column-major
 *              order to provide interchangability with the BLAS and LAPACK
 *              Fortran libraries.
 *
 * @tparam      eT An element type which represents a number that provides all common
 *              mathmatical operations.
 *
 * @since       0.0.1
 * 
 * @todo        Implement all other additional operators for 'lhs * rhs' outside
 *              of class definition
 *
 * @author      Denis-Michael Lux <denis.lux@icloud.com>
 * @date        13.05.15
 *
 * @ingroup     matrix
 */
template< typename eT >
class matrix< complex< eT > >
{
    size_t rows;                //!< Number of matrix rows
    size_t cols;                //!< Number of matrix columns
    
    size_t r_inj;               //!< Row injection index
    size_t c_inj;               //!< Column injection index
    
    complex< eT >* mem;         //!< Matrix memory
    
    // Friend classes to grant member access
    template< typename F > friend class matrix;
    template< typename F > friend class vector;
    
public:
    inline                               ~matrix();
    inline                                matrix();
    
    inline                                matrix(const size_t& m, const size_t& n);
    inline                                matrix(const size_t& mn);
    inline                                matrix(const size_t& m, const size_t& n, const complex< eT >& initial);
    
    inline                                matrix(const matrix< complex< eT > >& A);
    inline                                matrix(const matrix< eT >& A);
    inline                                matrix(matrix< complex< eT > >&& A);
    
    inline                                matrix(const vector< eT >& v);
    inline                                matrix(const vector< complex< eT > >& v);
    
    inline       matrix< complex< eT > >  operator+(const matrix< complex< eT > >& A);
    inline       matrix< complex< eT > >  operator+(const matrix< eT >& A);
    inline       matrix< complex< eT > >  operator-(const matrix< complex< eT > >& A);
    inline       matrix< complex< eT > >  operator-(const matrix< eT >& A);
    inline       matrix< complex< eT > >  operator*(const matrix< complex< eT > >& A);
    inline       matrix< complex< eT > >  operator*(const matrix< eT >& A);
    inline       matrix< complex< eT > >  operator%(const matrix< complex< eT > >& A);
    inline       matrix< complex< eT > >  operator%(const matrix< eT >& A);
    inline       matrix< complex< eT > >  operator/(const matrix< complex< eT > >& A);
    inline       matrix< complex< eT > >  operator/(const matrix< eT >& A);
    
    inline const matrix< complex< eT > >& operator=(const matrix< eT >& A);
    inline const matrix< complex< eT > >& operator=(const matrix< complex< eT > >& A);
    inline const matrix< complex< eT > >& operator=(matrix< complex< eT > >&& A);
    
    inline       bool                     operator==(const matrix< complex< eT > >& A);
    inline       bool                     operator!=(const matrix< complex< eT > >& A);
    
    inline const matrix< complex< eT > >& operator+=(const matrix< eT >& A);
    inline const matrix< complex< eT > >& operator+=(const matrix< complex< eT > >& A);
    inline const matrix< complex< eT > >& operator-=(const matrix< eT >& A);
    inline const matrix< complex< eT > >& operator-=(const matrix< complex< eT > >& A);
    inline const matrix< complex< eT > >& operator*=(const matrix< eT >& A);
    inline const matrix< complex< eT > >& operator*=(const matrix< complex< eT > >& A);
    inline const matrix< complex< eT > >& operator%=(const matrix< eT >& A);
    inline const matrix< complex< eT > >& operator%=(const matrix< complex< eT > >& A);
    inline const matrix< complex< eT > >& operator/=(const matrix< eT >& A);
    inline const matrix< complex< eT > >& operator/=(const matrix< complex< eT > >& A);
    
    inline const matrix< complex< eT > >& operator*=(const vector< eT >& v);
    inline const matrix< complex< eT > >& operator*=(const vector< complex< eT > >& v);
    
    inline       matrix< complex< eT > >  operator+();
    inline       matrix< complex< eT > >  operator-();
    
    inline       matrix< complex< eT > >  operator+(const eT& rhs);
    inline       matrix< complex< eT > >  operator+(const complex< eT >& rhs);
    inline       matrix< complex< eT > >  operator-(const eT& rhs);
    inline       matrix< complex< eT > >  operator-(const complex< eT >& rhs);
    inline       matrix< complex< eT > >  operator*(const eT& rhs);
    inline       matrix< complex< eT > >  operator*(const complex< eT >& rhs);
    inline       matrix< complex< eT > >  operator/(const eT& rhs);
    inline       matrix< complex< eT > >  operator/(const complex< eT >& rhs);
    inline       matrix< complex< eT > >  operator^(const unsigned int& exp);
    
    inline       vector< complex< eT > >  operator*(const vector< eT >& v);
    inline       vector< complex< eT > >  operator*(const vector< complex< eT > >& v);
    
    inline       matrix< complex< eT > >& operator+=(const eT& rhs);
    inline       matrix< complex< eT > >& operator+=(const complex< eT >& rhs);
    inline       matrix< complex< eT > >& operator-=(const eT& rhs);
    inline       matrix< complex< eT > >& operator-=(const complex< eT >& rhs);
    inline       matrix< complex< eT > >& operator*=(const eT& rhs);
    inline       matrix< complex< eT > >& operator*=(const complex< eT >& rhs);
    inline       matrix< complex< eT > >& operator/=(const eT& rhs);
    inline       matrix< complex< eT > >& operator/=(const complex< eT >& rhs);
    inline       matrix< complex< eT > >& operator^=(const unsigned int& exp);
    
    inline       matrix< complex< eT > >& operator<<(const eT& val);
    inline       matrix< complex< eT > >& operator<<(const complex< eT >& val);
    //inline matrix< eT >& operator<<(const tok& T);
    
    inline       complex< eT >&           operator()(const size_t& i, const size_t& j);
    inline const complex< eT >&           operator()(const size_t& i, const size_t& j) const;
    
    inline       void                     zeros();
    inline       void                     ones();
    inline       void                     eye();
    inline       void                     transpose();
    inline       void                     fill(const eT& value);
    inline       void                     fill(const complex< eT >& value);
    inline       void                     diag(const complex< eT >& value);
    inline       void                     diag(const vector< complex< eT > >& vec);
    
    inline constexpr size_t               n_rows() const;
    inline constexpr size_t               n_cols() const;
    inline       complex< eT >*           memptr();
    inline const complex< eT >*           memptr() const;
    
    inline const complex< double >        determinant();
};

template<typename eT> std::ostream& operator<<(std::ostream& o, const matrix< complex< eT > >& A);

UZLMATH_END

#endif
