//
//  vector_dec.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 05.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_vector_dec_hpp
#define UZLMathLib_vector_dec_hpp

UZLMATH_BEGIN

/*!
 * @brief       Collection of classes and functions for vectors for mathematical 
 *              purposes.
 * @defgroup    vector Vector
 * @{
 */
 
/*!
 * @brief       Vector class with special operations and methods.
 * @details     The vector class uses internally an dynamic allocated array to store 
 *              values. The memory is organized in as linear memory. It provides a 
 *              set of methods and operators to provide a useful vector class for 
 *              mathematical purposes.
 *
 *
 * @tparam      T An element type which represents a number that provides all common
 *              mathmatical operations.
 *
 * @since       0.1.1
 *
 * @author      Denis-Michael Lux <denis.lux@icloud.com>
 * @date        30.05.15
 */
template< typename T >
class vector
{
    size_t   inj;         //!< Index for injection
    
public:
    // ivars
    const size_t   size;  //!< size of vector
    const vec_type type;  //!< type of vector
    
    const T* mem;        //!< vector data
    
    // methods
    inline                               vector();
    inline                               vector(const size_t& s, const vec_type& t = vec_type::ROW);
    inline                               vector(const size_t& s, const T& initial, const vec_type& t = vec_type::ROW /*! default type */);
    inline                               vector(const vector< T >& vec);
    inline                               vector(const vector< T >& vec, const vec_type& t);
    inline                               vector(vector< T >&& vec);
    inline                              ~vector();
    
    inline       vector< T >             operator+(const vector< T >& v);
    inline       vector< T >             operator-(const vector< T >& v);
    inline       matrix< T >             operator*(const vector< T >& v);
    inline       vector< T >             operator/(const vector< T >& v);
    inline       vector< T >             operator%(const vector< T >& v);
    
    inline       vector< complex< T > >  operator+(const vector< complex< T > >& v);
    inline       vector< complex< T > >  operator-(const vector< complex< T > >& v);
    inline       matrix< complex< T > >  operator*(const vector< complex< T > >& v);
    inline       vector< complex< T > >  operator/(const vector< complex< T > >& v);
    inline       vector< complex< T > >  operator%(const vector< complex< T > >& v);
    
    inline       vector< T >             operator+(const T& s);
    inline       vector< T >             operator-(const T& s);
    inline       vector< T >             operator*(const T& s);
    inline       vector< T >             operator/(const T& s);
    
    inline       vector< complex< T > >  operator+(const complex< T >& s);
    inline       vector< complex< T > >  operator-(const complex< T >& s);
    inline       vector< complex< T > >  operator*(const complex< T >& s);
    inline       vector< complex< T > >  operator/(const complex< T >& s);
    
    inline       vector< T >             operator+();
    inline       vector< T >             operator-();
    
    inline       vector< T >             operator*(const matrix< T >& mat);
    inline       vector< complex< T > >  operator*(const matrix< complex< T > >& mat);
    
    inline       bool                    operator>(const vector< T >& v);
    inline       bool                    operator<(const vector< T >& v);
    
    inline const vector< T >&            operator=(const vector< T >& v);
    inline const vector< T >&            operator=(vector< T >&& v);
    
    inline const vector< T >&            operator+=(const vector< T >& v);
    inline const vector< T >&            operator-=(const vector< T >& v);
    inline const vector< T >&            operator/=(const vector< T >& v);
    inline const vector< T >&            operator%=(const vector< T >& v);
    inline const vector< T >&            operator*=(const vector< T >& v);
    
    inline const vector< T >&            operator+=(const T& s);
    inline const vector< T >&            operator-=(const T& s);
    inline const vector< T >&            operator*=(const T& s);
    inline const vector< T >&            operator/=(const T& s);
    
    inline       bool                    operator==(const vector< T >& v);
    inline       bool                    operator!=(const vector< T >& v);
    inline       bool                    operator>=(const vector< T >& v);
    inline       bool                    operator<=(const vector< T >& v);
    
    inline       T&                      operator[](const size_t& idx);
    inline const T&                      operator[](const size_t& idx) const;
    
    inline       T&                      operator()(const size_t& idx);
    inline const T&                      operator()(const size_t& idx) const;
    
    inline       void                    ones();
    inline       void                    zeros();
    inline       void                    transpose();
    inline       void                    fill(const T& s);
    
};

template< typename S > std::ostream& operator<<(std::ostream& o, const vector< S >& v);

/*!
 * @}
 */

UZLMATH_END

#endif
