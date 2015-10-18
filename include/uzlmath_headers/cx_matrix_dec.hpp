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
class
matrix< complex< T >, if_pod_type< T > >
{
    size_t r_inj;                   //!< Row injection index
    size_t c_inj;                   //!< Column injection index
    
public:
    // type declarations
    typedef T pod_type;             //!< The POD type of matrix elements
    
    // ivars
    const size_t rows;              //!< Number of matrix rows
    const size_t cols;              //!< Number of matrix columns
    
    const complex< pod_type >* mem; //!< Matrix memory
    
    // methods
    inline                                     ~matrix();
    inline                                      matrix();
    
    inline                                      matrix(const size_t& m, const size_t& n);
    inline                                      matrix(const size_t& mn);
    inline                                      matrix(const size_t& m, const size_t& n, const complex< pod_type >& initial);
    
    inline                                      matrix(const matrix< complex< pod_type > >& A);
    inline                                      matrix(const matrix< pod_type >& A);
    inline                                      matrix(matrix< complex< pod_type > >&& A);
    
    inline                                      matrix(const vector< pod_type >& v);
    inline                                      matrix(const vector< complex< pod_type > >& v);
    
    inline       matrix< complex< pod_type > >  operator+(const matrix< complex< pod_type > >& A);
    inline       matrix< complex< pod_type > >  operator+(const matrix< pod_type >& A);
    inline       matrix< complex< pod_type > >  operator-(const matrix< complex< pod_type > >& A);
    inline       matrix< complex< pod_type > >  operator-(const matrix< pod_type >& A);
    inline       matrix< complex< pod_type > >  operator*(const matrix< complex< pod_type > >& A);
    inline       matrix< complex< pod_type > >  operator*(const matrix< pod_type >& A);
    inline       matrix< complex< pod_type > >  operator%(const matrix< complex< pod_type > >& A);
    inline       matrix< complex< pod_type > >  operator%(const matrix< pod_type >& A);
    inline       matrix< complex< pod_type > >  operator/(const matrix< complex< pod_type > >& A);
    inline       matrix< complex< pod_type > >  operator/(const matrix< pod_type >& A);
    
    inline const matrix< complex< pod_type > >& operator=(const matrix< pod_type >& A);
    inline const matrix< complex< pod_type > >& operator=(const matrix< complex< pod_type > >& A);
    inline const matrix< complex< pod_type > >& operator=(matrix< complex< pod_type > >&& A);
    
    inline       bool                           operator==(const matrix< complex< pod_type > >& A);
    inline       bool                           operator!=(const matrix< complex< pod_type > >& A);
    
    inline const matrix< complex< pod_type > >& operator+=(const matrix< pod_type >& A);
    inline const matrix< complex< pod_type > >& operator+=(const matrix< complex< pod_type > >& A);
    inline const matrix< complex< pod_type > >& operator-=(const matrix< pod_type >& A);
    inline const matrix< complex< pod_type > >& operator-=(const matrix< complex< pod_type > >& A);
    inline const matrix< complex< pod_type > >& operator*=(const matrix< pod_type >& A);
    inline const matrix< complex< pod_type > >& operator*=(const matrix< complex< pod_type > >& A);
    inline const matrix< complex< pod_type > >& operator%=(const matrix< pod_type >& A);
    inline const matrix< complex< pod_type > >& operator%=(const matrix< complex< pod_type > >& A);
    inline const matrix< complex< pod_type > >& operator/=(const matrix< pod_type >& A);
    inline const matrix< complex< pod_type > >& operator/=(const matrix< complex< pod_type > >& A);
    
    inline const matrix< complex< pod_type > >& operator*=(const vector< pod_type >& v);
    inline const matrix< complex< pod_type > >& operator*=(const vector< complex< pod_type > >& v);
    
    inline       matrix< complex< pod_type > >  operator+();
    inline       matrix< complex< pod_type > >  operator-();
    
    inline       matrix< complex< pod_type > >  operator+(const pod_type& rhs);
    inline       matrix< complex< pod_type > >  operator+(const complex< pod_type >& rhs);
    inline       matrix< complex< pod_type > >  operator-(const pod_type& rhs);
    inline       matrix< complex< pod_type > >  operator-(const complex< pod_type >& rhs);
    inline       matrix< complex< pod_type > >  operator*(const pod_type& rhs);
    inline       matrix< complex< pod_type > >  operator*(const complex< pod_type >& rhs);
    inline       matrix< complex< pod_type > >  operator/(const pod_type& rhs);
    inline       matrix< complex< pod_type > >  operator/(const complex< pod_type >& rhs);
    inline       matrix< complex< pod_type > >  operator^(const unsigned int& exp);
    
    inline       vector< complex< pod_type > >  operator*(const vector< pod_type >& v);
    inline       vector< complex< pod_type > >  operator*(const vector< complex< pod_type > >& v);
    
    inline       matrix< complex< pod_type > >& operator+=(const pod_type& rhs);
    inline       matrix< complex< pod_type > >& operator+=(const complex< pod_type >& rhs);
    inline       matrix< complex< pod_type > >& operator-=(const pod_type& rhs);
    inline       matrix< complex< pod_type > >& operator-=(const complex< pod_type >& rhs);
    inline       matrix< complex< pod_type > >& operator*=(const pod_type& rhs);
    inline       matrix< complex< pod_type > >& operator*=(const complex< pod_type >& rhs);
    inline       matrix< complex< pod_type > >& operator/=(const pod_type& rhs);
    inline       matrix< complex< pod_type > >& operator/=(const complex< pod_type >& rhs);
    inline       matrix< complex< pod_type > >& operator^=(const unsigned int& exp);
    
    inline       matrix< complex< pod_type > >& operator<<(const pod_type& val);
    inline       matrix< complex< pod_type > >& operator<<(const complex< pod_type >& val);
    //inline matrix< T >& operator<<(const tok& T);
    
    inline       complex< pod_type >&           operator()(const size_t& i, const size_t& j);
    inline const complex< pod_type >&           operator()(const size_t& i, const size_t& j) const;
    
    inline       void                           zeros();
    inline       void                           ones();
    inline       void                           eye();
    inline       void                           transpose();
    inline       void                           fill(const pod_type& value);
    inline       void                           fill(const complex< pod_type >& value);
    inline       void                           diag(const complex< pod_type >& value);
    inline       void                           diag(const vector< complex< pod_type > >& vec);
    
    inline const complex< double >              determinant();
    
    // iterators
    class row_iterator
    {
        const pod_type* position;   //!< current position
        
    public:
        
    };
};

template<typename T> std::ostream& operator<<(std::ostream& o, const matrix< complex< T > >& A);

UZLMATH_END

#endif
