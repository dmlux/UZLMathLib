//
//  traits.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 11.01.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_traits_hpp
#define UZLMathLib_traits_hpp

UZLMATH_BEGIN

/*!
 * @brief       Collection of trait classes.
 * @details     Trait classes can be used to detect the type of template
 *              variables at compile time to optimize function efficiency.
 * @defgroup    traits Trait structs
 * @{
 */

// check for same type
template< typename T, typename U > struct same_type                     { static const bool value = false;  };
template< typename T >             struct same_type< T, T >             { static const bool value = true;   };

// check for different types
template< typename T, typename U > struct different_type                { static const bool value = true;   };
template< typename T >             struct different_type< T, T >        { static const bool value = false;  };

// check if type is complex number
template< typename T > struct is_complex                                { static const bool value = false;  };
template<>             struct is_complex< complex< short > >            { static const bool value = true;   };
template<>             struct is_complex< complex< int > >              { static const bool value = true;   };
template<>             struct is_complex< complex< long > >             { static const bool value = true;   };
template<>             struct is_complex< complex< long long > >        { static const bool value = true;   };
template<>             struct is_complex< complex< unsigned short > >   { static const bool value = true;   };
template<>             struct is_complex< complex< unsigned int > >     { static const bool value = true;   };
template<>             struct is_complex< complex< unsigned long > >    { static const bool value = true;   };
template<>             struct is_complex< complex< float > >            { static const bool value = true;   };
template<>             struct is_complex< complex< double > >           { static const bool value = true;   };
template<>             struct is_complex< complex< long double > >      { static const bool value = true;   };

// check if type is real
template< typename T > struct is_real_type                              { static const bool value = false;  };
template<>             struct is_real_type< float >                     { static const bool value = true;   };
template<>             struct is_real_type< double >                    { static const bool value = true;   };
template<>             struct is_real_type< long double >               { static const bool value = true;   };

// check if type is integral
template< typename T > struct is_integral_type                          { static const bool value = false;  };
template<>             struct is_integral_type< short >                 { static const bool value = true;   };
template<>             struct is_integral_type< int >                   { static const bool value = true;   };
template<>             struct is_integral_type< long >                  { static const bool value = true;   };
template<>             struct is_integral_type< long long >             { static const bool value = true;   };
template<>             struct is_integral_type< unsigned short >        { static const bool value = true;   };
template<>             struct is_integral_type< unsigned int >          { static const bool value = true;   };
template<>             struct is_integral_type< unsigned long >         { static const bool value = true;   };
template<>             struct is_integral_type< unsigned long long >    { static const bool value = true;   };

/*!
 * @}
 */

UZLMATH_END

#endif /* traits.hpp */
