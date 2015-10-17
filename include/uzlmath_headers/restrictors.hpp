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
template< typename T > struct uzl_void_real_only                                {                                               };
template<>             struct uzl_void_real_only< float >                       { typedef void result;                          };
template<>             struct uzl_void_real_only< double >                      { typedef void result;                          };
template<>             struct uzl_void_real_only< long double >                 { typedef void result;                          };

// Double with return type void
template< typename T > struct uzl_double_only                                   {                                               };
template<>             struct uzl_double_only< double >                         { typedef void result;                          };

// Real numbers for return type void
template< typename T > struct uzl_num_only                                      {                                               };
template<>             struct uzl_num_only< short >                             { typedef short result;                         };
template<>             struct uzl_num_only< int >                               { typedef int result;                           };
template<>             struct uzl_num_only< long >                              { typedef long result;                          };
template<>             struct uzl_num_only< long long >                         { typedef long long result;                     };
template<>             struct uzl_num_only< unsigned short >                    { typedef unsigned short result;                };
template<>             struct uzl_num_only< unsigned int >                      { typedef unsigned int result;                  };
template<>             struct uzl_num_only< unsigned long >                     { typedef unsigned long result;                 };
template<>             struct uzl_num_only< unsigned long long >                { typedef unsigned long long result;            };
template<>             struct uzl_num_only< float >                             { typedef float result;                         };
template<>             struct uzl_num_only< double >                            { typedef double result;                        };
template<>             struct uzl_num_only< long double >                       { typedef long double result;                   };

// Real numbers for return type void
template< typename T > struct uzl_void_num_only                                 {                                               };
template<>             struct uzl_void_num_only< short >                        { typedef void result;                          };
template<>             struct uzl_void_num_only< int >                          { typedef void result;                          };
template<>             struct uzl_void_num_only< long >                         { typedef void result;                          };
template<>             struct uzl_void_num_only< long long >                    { typedef void result;                          };
template<>             struct uzl_void_num_only< unsigned short >               { typedef void result;                          };
template<>             struct uzl_void_num_only< unsigned int >                 { typedef void result;                          };
template<>             struct uzl_void_num_only< unsigned long >                { typedef void result;                          };
template<>             struct uzl_void_num_only< unsigned long long >           { typedef void result;                          };
template<>             struct uzl_void_num_only< float >                        { typedef void result;                          };
template<>             struct uzl_void_num_only< double >                       { typedef void result;                          };
template<>             struct uzl_void_num_only< long double >                  { typedef void result;                          };

// Complex numbers for return type void
template< typename T > struct uzl_void_cx_num_only                              {                                               };
template<>             struct uzl_void_cx_num_only< complex< short > >          { typedef void result;                          };
template<>             struct uzl_void_cx_num_only< complex< int > >            { typedef void result;                          };
template<>             struct uzl_void_cx_num_only< complex< long > >           { typedef void result;                          };
template<>             struct uzl_void_cx_num_only< complex< long long > >      { typedef void result;                          };
template<>             struct uzl_void_cx_num_only< complex< unsigned short > > { typedef void result;                          };
template<>             struct uzl_void_cx_num_only< complex< unsigned int > >   { typedef void result;                          };
template<>             struct uzl_void_cx_num_only< complex< unsigned long > >  { typedef void result;                          };
template<>             struct uzl_void_cx_num_only< complex< float > >          { typedef void result;                          };
template<>             struct uzl_void_cx_num_only< complex< double > >         { typedef void result;                          };
template<>             struct uzl_void_cx_num_only< complex< long double > >    { typedef void result;                          };

// Real numbers for return type complex double vector
template< typename T > struct uzl_vec_cx_dbl_real_num_only                      {                                               };
template<>             struct uzl_vec_cx_dbl_real_num_only< short >             { typedef vector< complex< double > > result;   };
template<>             struct uzl_vec_cx_dbl_real_num_only< int >               { typedef vector< complex< double > > result;   };
template<>             struct uzl_vec_cx_dbl_real_num_only< long >              { typedef vector< complex< double > > result;   };
template<>             struct uzl_vec_cx_dbl_real_num_only< long long >         { typedef vector< complex< double > > result;   };
template<>             struct uzl_vec_cx_dbl_real_num_only< unsigned short >    { typedef vector< complex< double > > result;   };
template<>             struct uzl_vec_cx_dbl_real_num_only< unsigned int >      { typedef vector< complex< double > > result;   };
template<>             struct uzl_vec_cx_dbl_real_num_only< unsigned long >     { typedef vector< complex< double > > result;   };
template<>             struct uzl_vec_cx_dbl_real_num_only< float >             { typedef vector< complex< double > > result;   };
template<>             struct uzl_vec_cx_dbl_real_num_only< double >            { typedef vector< complex< double > > result;   };
template<>             struct uzl_vec_cx_dbl_real_num_only< long double >       { typedef vector< complex< double > > result;   };

// Integral random dist types
template< typename T > struct uzl_int_rand_dist_type                            { typedef int result;                           };
template<>             struct uzl_int_rand_dist_type< short >                   { typedef short result;                         };
template<>             struct uzl_int_rand_dist_type< int >                     { typedef int result;                           };
template<>             struct uzl_int_rand_dist_type< long >                    { typedef long result;                          };
template<>             struct uzl_int_rand_dist_type< long long >               { typedef long long result;                     };
template<>             struct uzl_int_rand_dist_type< unsigned short >          { typedef unsigned short result;                };
template<>             struct uzl_int_rand_dist_type< unsigned int >            { typedef unsigned int result;                  };
template<>             struct uzl_int_rand_dist_type< unsigned long >           { typedef unsigned long result;                 };
template<>             struct uzl_int_rand_dist_type< float >                   { typedef long result;                          };
template<>             struct uzl_int_rand_dist_type< double >                  { typedef long result;                          };
template<>             struct uzl_int_rand_dist_type< long double >             { typedef long long result;                     };

// Uniform real dist types
template< typename T > struct uzl_uniform_real_dist_only                                                {                       };
template<>             struct uzl_uniform_real_dist_only< uniform_real_distribution< float > >          { typedef void result;  };
template<>             struct uzl_uniform_real_dist_only< uniform_real_distribution< double > >         { typedef void result;  };
template<>             struct uzl_uniform_real_dist_only< uniform_real_distribution< long double > >    { typedef void result;  };

/*!
 * @}
 */

UZLMATH_END

#endif
