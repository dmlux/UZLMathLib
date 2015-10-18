//
//  restrictors.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 31.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_restrictors_hpp
#define UZLMathLib_restrictors_hpp

UZLMATH_BEGIN

/*!
 * @brief       Collection of restrictors to restrict element types for special functions
 * @defgroup    restrictors Restrictors
 * @{
 */

// Real values for return type void
template< typename T > struct void_real_only                                    {                                               };
template<>             struct void_real_only< float >                           { typedef void result;                          };
template<>             struct void_real_only< double >                          { typedef void result;                          };
template<>             struct void_real_only< long double >                     { typedef void result;                          };

template< typename T > using void_real_type = typename void_real_only< T >::result;

// Numbers only
template< typename T > struct numbers_only                                      {                                               };
template<>             struct numbers_only< short >                             { typedef short result;                         };
template<>             struct numbers_only< int >                               { typedef int result;                           };
template<>             struct numbers_only< long >                              { typedef long result;                          };
template<>             struct numbers_only< long long >                         { typedef long long result;                     };
template<>             struct numbers_only< unsigned short >                    { typedef unsigned short result;                };
template<>             struct numbers_only< unsigned int >                      { typedef unsigned int result;                  };
template<>             struct numbers_only< unsigned long >                     { typedef unsigned long result;                 };
template<>             struct numbers_only< unsigned long long >                { typedef unsigned long long result;            };
template<>             struct numbers_only< float >                             { typedef float result;                         };
template<>             struct numbers_only< double >                            { typedef double result;                        };
template<>             struct numbers_only< long double >                       { typedef long double result;                   };

template< typename T > using number_type = typename numbers_only< T >::result;

// Real numbers for return type void
template< typename T > struct void_numbers_only                                 {                                               };
template<>             struct void_numbers_only< short >                        { typedef void result;                          };
template<>             struct void_numbers_only< int >                          { typedef void result;                          };
template<>             struct void_numbers_only< long >                         { typedef void result;                          };
template<>             struct void_numbers_only< long long >                    { typedef void result;                          };
template<>             struct void_numbers_only< unsigned short >               { typedef void result;                          };
template<>             struct void_numbers_only< unsigned int >                 { typedef void result;                          };
template<>             struct void_numbers_only< unsigned long >                { typedef void result;                          };
template<>             struct void_numbers_only< unsigned long long >           { typedef void result;                          };
template<>             struct void_numbers_only< float >                        { typedef void result;                          };
template<>             struct void_numbers_only< double >                       { typedef void result;                          };
template<>             struct void_numbers_only< long double >                  { typedef void result;                          };

template< typename T > using void_number_type = typename void_numbers_only< T >::result;

// Complex numbers only
template< typename T > struct cx_numbers_only                                   {                                               };
template<>             struct cx_numbers_only< complex< short > >               { typedef complex< short > result;              };
template<>             struct cx_numbers_only< complex< int > >                 { typedef complex< int > result;                };
template<>             struct cx_numbers_only< complex< long > >                { typedef complex< long > result;               };
template<>             struct cx_numbers_only< complex< long long > >           { typedef complex< long long > result;          };
template<>             struct cx_numbers_only< complex< unsigned short > >      { typedef complex< unsigned short > result;     };
template<>             struct cx_numbers_only< complex< unsigned int > >        { typedef complex< unsigned int > result;       };
template<>             struct cx_numbers_only< complex< unsigned long > >       { typedef complex< unsigned long > result;      };
template<>             struct cx_numbers_only< complex< unsigned long long > >  { typedef complex< unsigned long long > result; };
template<>             struct cx_numbers_only< complex< float > >               { typedef complex< float > result;              };
template<>             struct cx_numbers_only< complex< double > >              { typedef complex< double > result;             };
template<>             struct cx_numbers_only< complex< long double > >         { typedef complex< long double > result;        };

template< typename T > using cx_number_type = typename cx_numbers_only< T >::result;

// Complex numbers for return type void
template< typename T > struct void_cx_numbers_only                              {                                               };
template<>             struct void_cx_numbers_only< complex< short > >          { typedef void result;                          };
template<>             struct void_cx_numbers_only< complex< int > >            { typedef void result;                          };
template<>             struct void_cx_numbers_only< complex< long > >           { typedef void result;                          };
template<>             struct void_cx_numbers_only< complex< long long > >      { typedef void result;                          };
template<>             struct void_cx_numbers_only< complex< unsigned short > > { typedef void result;                          };
template<>             struct void_cx_numbers_only< complex< unsigned int > >   { typedef void result;                          };
template<>             struct void_cx_numbers_only< complex< unsigned long > >  { typedef void result;                          };
template<>             struct void_cx_numbers_only< complex< float > >          { typedef void result;                          };
template<>             struct void_cx_numbers_only< complex< double > >         { typedef void result;                          };
template<>             struct void_cx_numbers_only< complex< long double > >    { typedef void result;                          };

template< typename T > using void_cx_number_type = typename void_cx_numbers_only< T >::result;

// Real numbers for return type complex double vector
template< typename T > struct cx_dblvec_numbers_only                            {                                               };
template<>             struct cx_dblvec_numbers_only< short >                   { typedef vector< complex< double > > result;   };
template<>             struct cx_dblvec_numbers_only< int >                     { typedef vector< complex< double > > result;   };
template<>             struct cx_dblvec_numbers_only< long >                    { typedef vector< complex< double > > result;   };
template<>             struct cx_dblvec_numbers_only< long long >               { typedef vector< complex< double > > result;   };
template<>             struct cx_dblvec_numbers_only< unsigned short >          { typedef vector< complex< double > > result;   };
template<>             struct cx_dblvec_numbers_only< unsigned int >            { typedef vector< complex< double > > result;   };
template<>             struct cx_dblvec_numbers_only< unsigned long >           { typedef vector< complex< double > > result;   };
template<>             struct cx_dblvec_numbers_only< float >                   { typedef vector< complex< double > > result;   };
template<>             struct cx_dblvec_numbers_only< double >                  { typedef vector< complex< double > > result;   };
template<>             struct cx_dblvec_numbers_only< long double >             { typedef vector< complex< double > > result;   };

template< typename T > using cx_dblvec_number_type = typename cx_dblvec_numbers_only< T >::result;

// Integral random dist types
template< typename T > struct int_rand_dist_type                                { typedef int result;                           };
template<>             struct int_rand_dist_type< short >                       { typedef short result;                         };
template<>             struct int_rand_dist_type< int >                         { typedef int result;                           };
template<>             struct int_rand_dist_type< long >                        { typedef long result;                          };
template<>             struct int_rand_dist_type< long long >                   { typedef long long result;                     };
template<>             struct int_rand_dist_type< unsigned short >              { typedef unsigned short result;                };
template<>             struct int_rand_dist_type< unsigned int >                { typedef unsigned int result;                  };
template<>             struct int_rand_dist_type< unsigned long >               { typedef unsigned long result;                 };
template<>             struct int_rand_dist_type< float >                       { typedef long result;                          };
template<>             struct int_rand_dist_type< double >                      { typedef long result;                          };
template<>             struct int_rand_dist_type< long double >                 { typedef long long result;                     };

template< typename T > using rand_int_dist_type = typename int_rand_dist_type< T >::result;

// Uniform real dist types
template< typename T > struct uniform_real_dist_only                                                {                           };
template<>             struct uniform_real_dist_only< uniform_real_distribution< float > >          { typedef void result;      };
template<>             struct uniform_real_dist_only< uniform_real_distribution< double > >         { typedef void result;      };
template<>             struct uniform_real_dist_only< uniform_real_distribution< long double > >    { typedef void result;      };

template< typename T > using uniform_real_dist_type = typename uniform_real_distribution< T >::result;

/*!
 * @}
 */

UZLMATH_END

#endif
