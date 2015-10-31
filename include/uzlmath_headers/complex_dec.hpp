//
//  comlex_dec.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 11.01.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_complex_dec_hpp
#define UZLMathLib_complex_dec_hpp

UZLMATH_BEGIN

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
 * @tparam  T An element type which represents a number that provides all common
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
template< typename T >
class complex
{
public:
    T re;   //!< The real part of the complex number
    T im;   //!< The imaginary part of the complex number
    
    inline                     ~complex();
    inline                      complex();
    
    inline                      complex(const T real_imag);
    inline                      complex(const T real, const T imag);
    
    inline                      complex(const complex< T >& c);
    
    inline const T              abs() const;
    inline const T              arg() const;
    inline const T              norm() const;
    
    inline       void           conj();
    inline       void           polar(const T& rho, const T& theta = 0);
    
    inline       complex< T >   operator+(const complex< T >& rhs);
    inline       complex< T >   operator-(const complex< T >& rhs);
    inline       complex< T >   operator*(const complex< T >& rhs);
    inline       complex< T >   operator/(const complex< T >& rhs);
    inline       complex< T >   operator^(const int rhs);
    
    inline const complex< T >&  operator+=(const complex< T >& rhs);
    inline const complex< T >&  operator-=(const complex< T >& rhs);
    inline const complex< T >&  operator*=(const complex< T >& rhs);
    inline const complex< T >&  operator/=(const complex< T >& rhs);
    inline const complex< T >&  operator^=(const int rhs);
    
    template< typename U > inline       complex< T >  operator+(const U& rhs);
    template< typename U > inline       complex< T >  operator-(const U& rhs);
    template< typename U > inline       complex< T >  operator*(const U& rhs);
    template< typename U > inline       complex< T >  operator/(const U& rhs);
    
    template< typename U > inline const complex< T >& operator+=(const U& rhs);
    template< typename U > inline const complex< T >& operator-=(const U& rhs);
    template< typename U > inline const complex< T >& operator*=(const U& rhs);
    template< typename U > inline const complex< T >& operator/=(const U& rhs);
};

/*- Additional operators -*/

//template< typename eT, typename T > complex< T >   operator*(const T& lhs, complex< T > rhs);
//template< typename eT, typename T > complex< T >   operator/(const T& lhs, complex< T > rhs);
//template< typename eT, typename T > complex< T >   operator+(const T& lhs, complex< T > rhs);
//template< typename eT, typename T > complex< T >   operator-(const T& lhs, complex< T > rhs);
//
//template< typename T >             complex< T >   operator+(int lhs, complex< T > rhs);

/*- Outstream operator -*/

template< typename T >             std::ostream&   operator<<(std::ostream& o, const complex< T >& c);

/*!
 * @}
 */

UZLMATH_END

#endif /* complex_dec.hpp */
