//
//  cx_matrix_dec.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 13.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_cx_matrix_dec_hpp
#define UZLMathLib_cx_matrix_dec_hpp

UZLMATH_BEGIN

/*!
 * @brief       Matrix specialization for complex data.
 * @details     A matrix that contains only complex values on each
 *              entry.  The matrix class uses internally an dynamic allocated array
 *              to store values. The memory is organized in  column-major
 *              order to provide interchangability with the BLAS and LAPACK
 *              Fortran libraries.
 *
 * @tparam      T An element type which represents a number that provides all common
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
template< typename T >
class matrix< complex< T > >
{
    size_t r_inj;               //!< Row injection index
    size_t c_inj;               //!< Column injection index
    
public:
    // ivars
    const size_t rows;          //!< Number of matrix rows
    const size_t cols;          //!< Number of matrix columns
    
    const complex< T >* mem;   //!< Matrix memory
    
    // methods
    inline                               ~matrix();
    inline                                matrix();
    
    inline                                matrix(const size_t& m, const size_t& n);
    inline                                matrix(const size_t& mn);
    inline                                matrix(const size_t& m, const size_t& n, const complex< T >& initial);
    
    inline                                matrix(const matrix< complex< T > >& A);
    inline                                matrix(const matrix< T >& A);
    inline                                matrix(matrix< complex< T > >&& A);
    
    inline                                matrix(const vector< T >& v);
    inline                                matrix(const vector< complex< T > >& v);
    
    inline       matrix< complex< T > >   operator+(const matrix< complex< T > >& A);
    inline       matrix< complex< T > >   operator+(const matrix< T >& A);
    inline       matrix< complex< T > >   operator-(const matrix< complex< T > >& A);
    inline       matrix< complex< T > >   operator-(const matrix< T >& A);
    inline       matrix< complex< T > >   operator*(const matrix< complex< T > >& A);
    inline       matrix< complex< T > >   operator*(const matrix< T >& A);
    inline       matrix< complex< T > >   operator%(const matrix< complex< T > >& A);
    inline       matrix< complex< T > >   operator%(const matrix< T >& A);
    inline       matrix< complex< T > >   operator/(const matrix< complex< T > >& A);
    inline       matrix< complex< T > >   operator/(const matrix< T >& A);
    
    inline const matrix< complex< T > >&  operator=(const matrix< T >& A);
    inline const matrix< complex< T > >&  operator=(const matrix< complex< T > >& A);
    inline const matrix< complex< T > >&  operator=(matrix< complex< T > >&& A);
    
    inline       bool                     operator==(const matrix< complex< T > >& A);
    inline       bool                     operator!=(const matrix< complex< T > >& A);
    
    inline const matrix< complex< T > >&  operator+=(const matrix< T >& A);
    inline const matrix< complex< T > >&  operator+=(const matrix< complex< T > >& A);
    inline const matrix< complex< T > >&  operator-=(const matrix< T >& A);
    inline const matrix< complex< T > >&  operator-=(const matrix< complex< T > >& A);
    inline const matrix< complex< T > >&  operator*=(const matrix< T >& A);
    inline const matrix< complex< T > >&  operator*=(const matrix< complex< T > >& A);
    inline const matrix< complex< T > >&  operator%=(const matrix< T >& A);
    inline const matrix< complex< T > >&  operator%=(const matrix< complex< T > >& A);
    inline const matrix< complex< T > >&  operator/=(const matrix< T >& A);
    inline const matrix< complex< T > >&  operator/=(const matrix< complex< T > >& A);
    
    inline const matrix< complex< T > >&  operator*=(const vector< T >& v);
    inline const matrix< complex< T > >&  operator*=(const vector< complex< T > >& v);
    
    inline       matrix< complex< T > >   operator+();
    inline       matrix< complex< T > >   operator-();
    
    inline       matrix< complex< T > >   operator+(const T& rhs);
    inline       matrix< complex< T > >   operator+(const complex< T >& rhs);
    inline       matrix< complex< T > >   operator-(const T& rhs);
    inline       matrix< complex< T > >   operator-(const complex< T >& rhs);
    inline       matrix< complex< T > >   operator*(const T& rhs);
    inline       matrix< complex< T > >   operator*(const complex< T >& rhs);
    inline       matrix< complex< T > >   operator/(const T& rhs);
    inline       matrix< complex< T > >   operator/(const complex< T >& rhs);
    inline       matrix< complex< T > >   operator^(const unsigned int& exp);
    
    inline       vector< complex< T > >   operator*(const vector< T >& v);
    inline       vector< complex< T > >   operator*(const vector< complex< T > >& v);
    
    inline       matrix< complex< T > >&  operator+=(const T& rhs);
    inline       matrix< complex< T > >&  operator+=(const complex< T >& rhs);
    inline       matrix< complex< T > >&  operator-=(const T& rhs);
    inline       matrix< complex< T > >&  operator-=(const complex< T >& rhs);
    inline       matrix< complex< T > >&  operator*=(const T& rhs);
    inline       matrix< complex< T > >&  operator*=(const complex< T >& rhs);
    inline       matrix< complex< T > >&  operator/=(const T& rhs);
    inline       matrix< complex< T > >&  operator/=(const complex< T >& rhs);
    inline       matrix< complex< T > >&  operator^=(const unsigned int& exp);
    
    inline       matrix< complex< T > >&  operator<<(const T& val);
    inline       matrix< complex< T > >&  operator<<(const complex< T >& val);
    //inline matrix< T >& operator<<(const tok& T);
    
    inline       complex< T >&            operator()(const size_t& i, const size_t& j);
    inline const complex< T >&            operator()(const size_t& i, const size_t& j) const;
    
    inline       void                     zeros();
    inline       void                     ones();
    inline       void                     eye();
    inline       void                     transpose();
    inline       void                     fill(const T& value);
    inline       void                     fill(const complex< T >& value);
    inline       void                     diag(const complex< T >& value);
    inline       void                     diag(const vector< complex< T > >& vec);
    
    inline const complex< double >        determinant();
};

template<typename T> std::ostream& operator<<(std::ostream& o, const matrix< complex< T > >& A);

UZLMATH_END

#endif
