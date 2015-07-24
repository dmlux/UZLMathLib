//
//  compiler.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 24.07.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_compiler_hpp
#define uzlmath_compiler_hpp

// Define the default number of used threads for multithreaded parts
#ifdef _OPENMP

    // include OpenMP headers for the OpenMP API
    #include <omp.h>

    // Number of default used threads is the max
    // number of possible threads
    #define UZL_MAX_THREADS omp_get_max_threads()
#else

    // setting the number of default used threads for
    // multithreaded parts to 1
    #define UZL_MAX_THREADS 1
#endif

/*- Namespace macros -*/
#define UZLMATH_BEGIN           namespace uzlmath {
#define UZLMATH_END             }

#define UZLMATH_NAMESPACE(a)    namespace uzlmath { namespace a {
#define UZLMATH_NAMESPACE_END   }}

/*- Debugging messages -*/


#endif
