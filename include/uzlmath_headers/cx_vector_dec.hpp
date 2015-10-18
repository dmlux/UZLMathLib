//
//  cx_vector_dec.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 02.06.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_cx_vector_dec_hpp
#define UZLMathLib_cx_vector_dec_hpp

UZLMATH_BEGIN

/*!
 * @ingroup    vector
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
vector< complex< T >, typename if_true< is_num_type< T >::value >::type >
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
    const complex< T >* mem;            //!< vector data
    
    inline                                vector();
    inline                                vector(const size_t& s, const vec_type& type = vec_type::ROW);
    inline                                vector(const size_t& s, const T& initial, const vec_type& type = vec_type::ROW);
    inline                                vector(const size_t& s, const complex< T >& initial, const vec_type& type = vec_type::ROW);
    inline                                vector(const vector< T >& vec);
    inline                                vector(const vector< complex< T > >& vec);
    inline                                vector(const vector< T >& vec, const vec_type& type);
    inline                                vector(const vector< complex< T > >& vec, const vec_type& type);
    inline                                vector(vector< complex< T > >&& vec);
    inline                               ~vector();
    
    inline       vector< complex< T > >   operator+(const vector< T >& v);
    inline       vector< complex< T > >   operator-(const vector< T >& v);
    inline       matrix< complex< T > >   operator*(const vector< T >& v);
    inline       vector< complex< T > >   operator/(const vector< T >& v);
    inline       vector< complex< T > >   operator%(const vector< T >& v);
    
    inline       vector< complex< T > >   operator+(const vector< complex< T > >& v);
    inline       vector< complex< T > >   operator-(const vector< complex< T > >& v);
    inline       matrix< complex< T > >   operator*(const vector< complex< T > >& v);
    inline       vector< complex< T > >   operator/(const vector< complex< T > >& v);
    inline       vector< complex< T > >   operator%(const vector< complex< T > >& v);
    
    inline       vector< complex< T > >   operator+(const T& s);
    inline       vector< complex< T > >   operator-(const T& s);
    inline       vector< complex< T > >   operator*(const T& s);
    inline       vector< complex< T > >   operator/(const T& s);
    
    inline       vector< complex< T > >   operator+(const complex< T >& s);
    inline       vector< complex< T > >   operator-(const complex< T >& s);
    inline       vector< complex< T > >   operator*(const complex< T >& s);
    inline       vector< complex< T > >   operator/(const complex< T >& s);
    
    inline       vector< complex< T > >   operator+();
    inline       vector< complex< T > >   operator-();
    
    inline       vector< complex< T > >   operator*(const matrix< T >& mat);
    inline       vector< complex< T > >   operator*(const matrix< complex< T > >& mat);
    
    inline const vector< complex< T > >&  operator=(const vector< T >& v);
    inline const vector< complex< T > >&  operator=(const vector< complex< T > >& v);
    inline const vector< complex< T > >&  operator=(vector< complex< T > >&& v);
    
    inline const vector< complex< T > >&  operator+=(const vector< T >& v);
    inline const vector< complex< T > >&  operator+=(const vector< complex< T > >& v);
    inline const vector< complex< T > >&  operator-=(const vector< T >& v);
    inline const vector< complex< T > >&  operator-=(const vector< complex< T > >& v);
    inline const vector< complex< T > >&  operator/=(const vector< T >& v);
    inline const vector< complex< T > >&  operator/=(const vector< complex< T > >& v);
    inline const vector< complex< T > >&  operator%=(const vector< T >& v);
    inline const vector< complex< T > >&  operator%=(const vector< complex< T > >& v);
    inline const vector< complex< T > >&  operator*=(const vector< T >& v);
    inline const vector< complex< T > >&  operator*=(const vector< complex< T > >& v);
    
    inline const vector< complex< T > >&  operator+=(const T& s);
    inline const vector< complex< T > >&  operator-=(const T& s);
    inline const vector< complex< T > >&  operator*=(const T& s);
    inline const vector< complex< T > >&  operator/=(const T& s);
    
    inline       bool                     operator==(const vector< T >& v);
    inline       bool                     operator!=(const vector< T >& v);
    inline       bool                     operator==(const vector< complex< T > >& v);
    inline       bool                     operator!=(const vector< complex< T > >& v);
    
    inline       complex< T >&            operator[](const size_t& idx);
    inline constexpr complex< T >&        operator[](const size_t& idx) const;
    
    inline       complex< T >&            operator()(const size_t& idx);
    inline constexpr complex< T >&        operator()(const size_t& idx) const;
    
    inline       void                     ones();
    inline       void                     zeros();
    inline       void                     transpose();
    inline       void                     fill(const T& s);
    
};

template< typename S > std::ostream& operator<<(std::ostream& o, const vector< complex< S > >& v);

/*!
 * @}
 */

UZLMATH_END

#endif
