//
//  comlex_dec.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 11.01.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_complex_dec_hpp
#define uzlmath_complex_dec_hpp

/*!
 * @brief       Collection of classes and functions for complex numbers for
 *              mathematical purposes.
 * @details     Contains a set of classes and functions that are useful for 
 *              dealing with complex numbers.
 * @defgroup    complex Complex
 * @{
 */

/*! 
 * @brief   A C++ implementation of a complex number
 * @details A complex number that contains all basic operations for calculating 
 *          with this type of numbers. All necessary operators are overloaded to 
 *          provide easier readability in algorihtms. The complex number uses a 
 *          template parameter to store the real and imaginary part in a member 
 *          variable of a specific type. Any POD type of C or C++ can be used to 
 *          store those parts of the complex number. Additionally to the build in 
 *          operators the calculation between PODs and a complex number are 
 *          implemented to.
 *
 * @tparam  eT An element type which represents a number that provides all common
 *          mathmatical operations.
 *
 * @since   0.0.1
 *
 * @todo    Add additional operators for different datatypes.
 * @todo    Write conversion operator overloads.
 *
 * @author  Denis-Michael Lux <denis.lux@icloud.com>
 * @date    11.01.15
 */
template< typename eT >
class complex
{
public:
    eT re;   //!< The real part of the complex number
    eT im;   //!< The imaginary part of the complex number
    
    inline                     ~complex();
    inline                      complex();
    
    inline                      complex(const eT real_imag);
    inline                      complex(const eT real, const eT imag);
    
    inline                      complex(const complex< eT >& c);
    
    inline const eT             abs() const;
    inline const eT             arg() const;
    inline const eT             norm() const;
    
    inline       void           conj();
    inline       void           polar(const eT& rho, const eT& theta = 0);
    
    inline       complex< eT >  operator+(const complex< eT >& rhs);
    inline       complex< eT >  operator-(const complex< eT >& rhs);
    inline       complex< eT >  operator*(const complex< eT >& rhs);
    inline       complex< eT >  operator/(const complex< eT >& rhs);
    inline       complex< eT >  operator^(const int rhs);
    
    inline const complex< eT >& operator+=(const complex< eT >& rhs);
    inline const complex< eT >& operator-=(const complex< eT >& rhs);
    inline const complex< eT >& operator*=(const complex< eT >& rhs);
    inline const complex< eT >& operator/=(const complex< eT >& rhs);
    inline const complex< eT >& operator^=(const int rhs);
    
    template< typename T > inline       complex< eT >  operator+(const T& rhs);
    template< typename T > inline       complex< eT >  operator-(const T& rhs);
    template< typename T > inline       complex< eT >  operator*(const T& rhs);
    template< typename T > inline       complex< eT >  operator/(const T& rhs);
    
    template< typename T > inline const complex< eT >& operator+=(const T& rhs);
    template< typename T > inline const complex< eT >& operator-=(const T& rhs);
    template< typename T > inline const complex< eT >& operator*=(const T& rhs);
    template< typename T > inline const complex< eT >& operator/=(const T& rhs);
};

/*- Additional operators -*/

template< typename eT, typename T > complex< eT >   operator*(const T& lhs, complex< eT > rhs);
template< typename eT, typename T > complex< eT >   operator/(const T& lhs, complex< eT > rhs);
template< typename eT, typename T > complex< eT >   operator+(const T& lhs, complex< eT > rhs);
template< typename eT, typename T > complex< eT >   operator-(const T& lhs, complex< eT > rhs);

template< typename eT >             complex< eT >   operator+(int lhs, complex< eT > rhs);

/*- Outstream operator -*/

template< typename eT >             std::ostream&   operator<<(std::ostream& o, const complex< eT >& c);

/*!
 * @}
 */

#endif
