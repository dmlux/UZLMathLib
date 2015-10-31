//
//  typedefs.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 25.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_typedefs_hpp
#define UZLMathLib_typedefs_hpp

UZLMATH_BEGIN

// matrix typedefs
typedef matrix< short >               s_mat;
typedef matrix< int >                 i_mat;
typedef matrix< long >                l_mat;
typedef matrix< long long >           ll_mat;
typedef matrix< unsigned short >      us_mat;
typedef matrix< unsigned int >        ui_mat;
typedef matrix< unsigned long >       ul_mat;
typedef matrix< unsigned long long >  ull_mat;
typedef matrix< float >               f_mat;
typedef matrix< double >              d_mat;
typedef matrix< long double >         ld_mat;
typedef matrix< complex< double > >   zx_mat;
typedef matrix< complex< float > >    cx_mat;

// vector typedefs
typedef vector< short >               s_vec;
typedef vector< int >                 i_vec;
typedef vector< long >                l_vec;
typedef vector< long long >           ll_vec;
typedef vector< unsigned short >      us_vec;
typedef vector< unsigned int >        ui_vec;
typedef vector< unsigned long >       ul_vec;
typedef vector< unsigned long long >  ull_vec;
typedef vector< float >               f_vec;
typedef vector< double >              d_vec;
typedef vector< long double >         ld_vec;
typedef vector< complex< double > >   zx_vec;
typedef vector< complex< float > >    cx_vec;

// grid3D typedefs
typedef grid3D< short >               s_grid3D;
typedef grid3D< int >                 i_grid3D;
typedef grid3D< long >                l_grid3D;
typedef grid3D< long long >           ll_grid3D;
typedef grid3D< unsigned short >      us_grid3D;
typedef grid3D< unsigned int >        ui_grid3D;
typedef grid3D< unsigned long >       ul_grid3D;
typedef grid3D< unsigned long long >  ull_grid3D;
typedef grid3D< float >               f_grid3D;
typedef grid3D< double >              d_grid3D;
typedef grid3D< long double >         ld_grid3D;
typedef grid3D< complex< double > >   zx_grid3D;
typedef grid3D< complex< float > >    cx_grid3D;

// factorial typedef
typedef factorial                     fac;

// complex typedefs
typedef complex< double >             complex_double;
typedef complex< double >             cx_double;
typedef complex< double >             zx;
typedef complex< float >              complex_float;
typedef complex< float >              cx_float;
typedef complex< float >              cx;

UZLMATH_END

#endif /* typedefs.hpp */
