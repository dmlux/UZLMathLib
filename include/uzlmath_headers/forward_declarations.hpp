//
//  forward_declarations.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 25.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_forward_declarations_hpp
#define UZLMathLib_forward_declarations_hpp

UZLMATH_BEGIN

// check if type is num type
template< typename T > struct is_num_type                               { static const bool value = false;  };
template<>             struct is_num_type< short >                      { static const bool value = true;   };
template<>             struct is_num_type< int >                        { static const bool value = true;   };
template<>             struct is_num_type< long >                       { static const bool value = true;   };
template<>             struct is_num_type< long long >                  { static const bool value = true;   };
template<>             struct is_num_type< unsigned short >             { static const bool value = true;   };
template<>             struct is_num_type< unsigned int >               { static const bool value = true;   };
template<>             struct is_num_type< unsigned long >              { static const bool value = true;   };
template<>             struct is_num_type< unsigned long long >         { static const bool value = true;   };
template<>             struct is_num_type< float >                      { static const bool value = true;   };
template<>             struct is_num_type< double >                     { static const bool value = true;   };
template<>             struct is_num_type< long double >                { static const bool value = true;   };

// check if template is true
template< bool T, typename U = void > struct if_true                    {                                   };
template< typename U >                struct if_true< true, U >         { typedef U type;                   };

// forward declarations
template< typename derived >                struct base;
template< typename T1, typename T2 >        class  glue;
template< typename T >                      class  complex;

// matrix
template< typename, typename = void >       class  matrix;
template< typename T >                      class  matrix< T,            typename if_true< is_num_type< T >::value >::type >;
template< typename T >                      class  matrix< complex< T >, typename if_true< is_num_type< T >::value >::type >;

// vector
template< typename, typename = void >       class  vector;
template< typename T >                      class  vector< T,            typename if_true< is_num_type< T >::value >::type >;
template< typename T >                      class  vector< complex< T >, typename if_true< is_num_type< T >::value >::type >;

template< typename T >                      struct grid3D;

template< typename T >                      class  memory;

template< typename T >                      class  array;

                                            class  access;

                                            class  factorial;
                                            class  stopwatch;

                                            struct DSOFTFourierCoefficients;

template< typename pod_type, typename derived > struct randctx;
template< typename T, typename A = void >       struct uniform_int_distribution;
template< typename T, typename A = void >       struct uniform_real_distribution;
template< typename T, typename A = void >       struct normal_distribution;

/*!
 * @brief       Specifing the type a vector can have
 * @details     Indicator for the type of the current vector.
 *              A vector can either be a column vector or a
 *              row vector.
 *
 * @ingroup     vector
 */
enum class vec_type
{
    ROW,        //!< Represents a column vector of dimension \f$M\times 1\f$ where \f$M\in\mathbb{N}^+\f$
    COLUMN      //!< Represents a row vector of dimension \f$1\times M\f$ where \f$M\in\mathbb{N}^+\f$
};

/*!
 * @brief       Injection operator tokens for matrices
 * @details     This enum class contains tokens to manipulate 
 *              the injection operator for classes.
 */
enum class mat
{
    NEXT_ROW    //!< Represents a token indicating the next row for injection operators
};

/*!
 * @brief       A collection of usable random engines
 */
enum class random_engine
{
    DEFAULT,            //!< Mersenne Twister 64 bit generator
    MINSTD_RAND,        //!< Minimal Standard generator
    MINSTD_RAND0,       //!< Minimal Standard 0 generator
    MERSENNE_TWISTER,   //!< Mersenne Twister generator
    MERSENNE_TWISTER64, //!< Mersenne Twister 64 bit generator
    RANLUX24_BASE,      //!< Ranlux 24 Base generator
    RANLUX48_BASE,      //!< Ranlux 48 Base generator
    RANLUX24,           //!< Ranlux 24 generator
    RANLUX48,           //!< Ranlux 48 generator
    KNUTH_B             //!< Knuth B generator
};

UZLMATH_END

#endif
