//
//  traits.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 11.01.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_traits_hpp
#define uzlmath_traits_hpp

/*!
 * @brief       Collection of trait classes.
 * @details     Trait classes can be used to detect the type of template
 *              variables at compile time to optimize function efficiency.
 * @defgroup    traits Trait structs
 * @{
 */

/*! @brief Trait that contains false value for each type that is not float. */
template< typename T > struct is_float                    { static const bool value = false; /*!< "is type of" flag */ };
/*! @brief Trait that contains true value for type float. */
template<>             struct is_float< float >           { static const bool value = true;  /*!< "is type of" flag */ };

/*! @brief Trait that contains false value for each type that is not double. */
template< typename T > struct is_double                   { static const bool value = false; /*!< "is type of" flag */ };
/*! @brief Trait that contains true value for type double. */
template<>             struct is_double< double >         { static const bool value = true;  /*!< "is type of" flag */ };

/*! @brief Trait that contains false value for each type that is not long double. */
template< typename T > struct is_ldouble                  { static const bool value = false; /*!< "is type of" flag */ };
/*! @brief Trait that contains true value for type long double. */
template<>             struct is_ldouble< long double >   { static const bool value = true;  /*!< "is type of" flag */ };

/*! @brief Trait that contains false value for each type that is not short. */
template< typename T > struct is_short                    { static const bool value = false; /*!< "is type of" flag */ };
/*! @brief Trait that contains true value for type short. */
template<>             struct is_short< short >           { static const bool value = true;  /*!< "is type of" flag */ };

/*! @brief Trait that contains false value for each type that is not int. */
template< typename T > struct is_int                      { static const bool value = false; /*!< "is type of" flag */ };
/*! @brief Trait that contains true value for type int. */
template<>             struct is_int< int >               { static const bool value = true;  /*!< "is type of" flag */ };

/*! @brief Trait that contains false value for each type that is not long. */
template< typename T > struct is_long                     { static const bool value = false; /*!< "is type of" flag */ };
/*! @brief Trait that contains true value for type long. */
template<>             struct is_long< long >             { static const bool value = true;  /*!< "is type of" flag */ };

/*! @brief Trait that contains false value for each type that is not unsigned short. */
template< typename T > struct is_ushort                   { static const bool value = false; /*!< "is type of" flag */ };
/*! @brief Trait that contains true value for type unsigned short. */
template<>             struct is_ushort< unsigned short > { static const bool value = true;  /*!< "is type of" flag */ };

/*! @brief Trait that contains false value for each type that is not unsigned int. */
template< typename T > struct is_uint                     { static const bool value = false; /*!< "is type of" flag */ };
/*! @brief Trait that contains true value for type unsigned int. */
template<>             struct is_uint< unsigned int >     { static const bool value = true;  /*!< "is type of" flag */ };

/*! @brief Trait that contains false value for each type that is not unsigned long. */
template< typename T > struct is_ulong                    { static const bool value = false; /*!< "is type of" flag */ };
/*! @brief Trait that contains true value for type unsigned long. */
template<>             struct is_ulong< unsigned long >   { static const bool value = true;  /*!< "is type of" flag */ };

/*! @brief Trait that contains false value for each type that is not long long. */
template< typename T > struct is_llong                    { static const bool value = false; /*!< "is type of" flag */ };
/*! @brief Trait that contains true value for type long long. */
template<>             struct is_llong< long long >       { static const bool value = true;  /*!< "is type of" flag */ };

/*! @brief Trait that contains false value for each type that is not complex. */
template< typename T > struct is_complex                  { static const bool value = false; /*!< "is type of" flag */ };
/*! @brief Trait that contains true value for type complex. */
template< typename T > struct is_complex< complex< T > >  { static const bool value = true;  /*!< "is type of" flag */ };

/*!
 * @}
 */

#endif
