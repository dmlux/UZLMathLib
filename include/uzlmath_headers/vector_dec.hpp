//
//  vector_dec.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 05.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_vector_dec_hpp
#define uzlmath_vector_dec_hpp

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
 * @tparam      eT An element type which represents a number that provides all common
 *              mathmatical operations.
 *
 * @since       0.1.1
 *
 * @author      Denis-Michael Lux <denis.lux@icloud.com>
 * @date        30.05.15
 */
template< typename eT >
class vector
{
private:
    
    size_t   size;  //!< size of vector
    vec_type type;  //!< type of vector
    
    size_t   inj;   //!< Index for injection
    
    eT* mem;        //!< vector data
    
    // Friend classes to grant member access
    template< typename F > friend class vector;
    template< typename F > friend class matrix;
    
public:
    inline                               vector();
    inline                               vector(const size_t& s, const vec_type& t = vec_type::ROW);
    inline                               vector(const size_t& s, const eT& initial, const vec_type& t = vec_type::ROW /*! default type */);
    inline                               vector(const vector< eT >& vec);
    inline                               vector(const vector< eT >& vec, const vec_type& t);
    inline                               vector(vector< eT >&& vec);
    inline                              ~vector();
    
    inline       vector< eT >            operator+(const vector< eT >& v);
    inline       vector< eT >            operator-(const vector< eT >& v);
    inline       matrix< eT >            operator*(const vector< eT >& v);
    inline       vector< eT >            operator/(const vector< eT >& v);
    inline       vector< eT >            operator%(const vector< eT >& v);
    
    inline       vector< complex< eT > > operator+(const vector< complex< eT > >& v);
    inline       vector< complex< eT > > operator-(const vector< complex< eT > >& v);
    inline       matrix< complex< eT > > operator*(const vector< complex< eT > >& v);
    inline       vector< complex< eT > > operator/(const vector< complex< eT > >& v);
    inline       vector< complex< eT > > operator%(const vector< complex< eT > >& v);
    
    inline       vector< eT >            operator+(const eT& s);
    inline       vector< eT >            operator-(const eT& s);
    inline       vector< eT >            operator*(const eT& s);
    inline       vector< eT >            operator/(const eT& s);
    
    inline       vector< complex< eT > > operator+(const complex< eT >& s);
    inline       vector< complex< eT > > operator-(const complex< eT >& s);
    inline       vector< complex< eT > > operator*(const complex< eT >& s);
    inline       vector< complex< eT > > operator/(const complex< eT >& s);
    
    inline       vector< eT >            operator+();
    inline       vector< eT >            operator-();
    
    inline       vector< eT >            operator*(const matrix< eT >& mat);
    inline       vector< complex< eT > > operator*(const matrix< complex< eT > >& mat);
    
    inline       bool                    operator>(const vector< eT >& v);
    inline       bool                    operator<(const vector< eT >& v);
    
    inline const vector< eT >&           operator=(const vector< eT >& v);
    inline const vector< eT >&           operator=(vector< eT >&& v);
    
    inline const vector< eT >&           operator+=(const vector< eT >& v);
    inline const vector< eT >&           operator-=(const vector< eT >& v);
    inline const vector< eT >&           operator/=(const vector< eT >& v);
    inline const vector< eT >&           operator%=(const vector< eT >& v);
    inline const vector< eT >&           operator*=(const vector< eT >& v);
    
    inline const vector< eT >&           operator+=(const eT& s);
    inline const vector< eT >&           operator-=(const eT& s);
    inline const vector< eT >&           operator*=(const eT& s);
    inline const vector< eT >&           operator/=(const eT& s);
    
    inline       bool                    operator==(const vector< eT >& v);
    inline       bool                    operator!=(const vector< eT >& v);
    inline       bool                    operator>=(const vector< eT >& v);
    inline       bool                    operator<=(const vector< eT >& v);
    
    inline       eT&                     operator[](const size_t& idx);
    inline constexpr eT&                 operator[](const size_t& idx) const;
    
    inline       eT&                     operator()(const size_t& idx);
    inline constexpr eT&                 operator()(const size_t& idx) const;
    
    inline       void                    ones();
    inline       void                    zeros();
    inline       void                    transpose();
    inline       void                    fill(const eT& s);
    
    inline       eT*                     memptr();
    inline const eT*                     memptr() const;
    
    inline constexpr size_t              n_elements() const;
    inline       vec_type                vecType() const;
    
};

template< typename S > std::ostream& operator<<(std::ostream& o, const vector< S >& v);

/*!
 * @}
 */

#endif
