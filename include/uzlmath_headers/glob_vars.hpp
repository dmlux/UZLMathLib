//
//  glob_vars.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 06.06.15.
//
//

#ifndef uzlmath_glob_vars_hpp
#define uzlmath_glob_vars_hpp

/*!
 * @brief           The global seed variable for the UZLMath library
 * @details         The seed is used for alle random functionality that is 
 *                  used in the UZLMath lib. Every time this variable is
 *                  used it gets updated with a new seed.
 */
static unsigned int uzlmath_seed = static_cast< unsigned int >(time(NULL));

#endif
