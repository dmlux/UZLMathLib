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
class
vector< T, if_pod_type< T > >
{
    size_t   inj;                       //!< Index for injection
    
public:
    // type declarations
    typedef T pod_type;                 //!< The POD type of matrix elements
    
    // iterators
    typedef const pod_type* iterator;   //!< Vector iterator
    
    // ivars
    const size_t   size;                //!< size of vector
    const vec_type type;                //!< type of vector
    const pod_type* mem;                //!< vector data
    
    // methods
    inline                                      vector();
    inline                                      vector(const size_t& s, const vec_type& t = vec_type::ROW);
    inline                                      vector(const size_t& s, const pod_type& initial, const vec_type& t = vec_type::ROW);
    inline                                      vector(const vector< pod_type >& vec);
    inline                                      vector(const vector< pod_type >& vec, const vec_type& t);
    inline                                      vector(vector< pod_type >&& vec);
    inline                                     ~vector();
    
    inline       vector< pod_type >             operator+(const vector< pod_type >& v);
    inline       vector< pod_type >             operator-(const vector< pod_type >& v);
    inline       matrix< pod_type >             operator*(const vector< pod_type >& v);
    inline       vector< pod_type >             operator/(const vector< pod_type >& v);
    inline       vector< pod_type >             operator%(const vector< pod_type >& v);
    
    inline       vector< complex< pod_type > >  operator+(const vector< complex< pod_type > >& v);
    inline       vector< complex< pod_type > >  operator-(const vector< complex< pod_type > >& v);
    inline       matrix< complex< pod_type > >  operator*(const vector< complex< pod_type > >& v);
    inline       vector< complex< pod_type > >  operator/(const vector< complex< pod_type > >& v);
    inline       vector< complex< pod_type > >  operator%(const vector< complex< pod_type > >& v);
    
    inline       vector< pod_type >             operator+(const pod_type& s);
    inline       vector< pod_type >             operator-(const pod_type& s);
    inline       vector< pod_type >             operator*(const pod_type& s);
    inline       vector< pod_type >             operator/(const pod_type& s);
    
    inline       vector< complex< pod_type > >  operator+(const complex< pod_type >& s);
    inline       vector< complex< pod_type > >  operator-(const complex< pod_type >& s);
    inline       vector< complex< pod_type > >  operator*(const complex< pod_type >& s);
    inline       vector< complex< pod_type > >  operator/(const complex< pod_type >& s);
    
    inline       vector< pod_type >             operator+();
    inline       vector< pod_type >             operator-();
    
    inline       vector< pod_type >             operator*(const matrix< pod_type >& mat);
    inline       vector< complex< pod_type > >  operator*(const matrix< complex< pod_type > >& mat);
    
    inline       bool                           operator>(const vector< pod_type >& v);
    inline       bool                           operator<(const vector< pod_type >& v);
    
    inline const vector< pod_type >&            operator=(const vector< pod_type >& v);
    inline const vector< pod_type >&            operator=(vector< pod_type >&& v);
    
    inline const vector< pod_type >&            operator+=(const vector< pod_type >& v);
    inline const vector< pod_type >&            operator-=(const vector< pod_type >& v);
    inline const vector< pod_type >&            operator/=(const vector< pod_type >& v);
    inline const vector< pod_type >&            operator%=(const vector< pod_type >& v);
    inline const vector< pod_type >&            operator*=(const vector< pod_type >& v);
    
    inline const vector< pod_type >&            operator+=(const pod_type& s);
    inline const vector< pod_type >&            operator-=(const pod_type& s);
    inline const vector< pod_type >&            operator*=(const pod_type& s);
    inline const vector< pod_type >&            operator/=(const pod_type& s);
    
    inline       bool                           operator==(const vector< pod_type >& v);
    inline       bool                           operator!=(const vector< pod_type >& v);
    inline       bool                           operator>=(const vector< pod_type >& v);
    inline       bool                           operator<=(const vector< pod_type >& v);
    
    inline       pod_type&                      operator[](const size_t& idx);
    inline const pod_type&                      operator[](const size_t& idx) const;
    
    inline       pod_type&                      operator()(const size_t& idx);
    inline const pod_type&                      operator()(const size_t& idx) const;
    
    inline       void                           ones();
    inline       void                           zeros();
    inline       void                           transpose();
    inline       void                           fill(const pod_type& s);
    
};

template< typename S > std::ostream& operator<<(std::ostream& o, const vector< S >& v);

/*!
 * @}
 */

UZLMATH_END

#endif
