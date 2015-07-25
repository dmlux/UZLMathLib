//
//  stop_watch_dec.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 13.01.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_stop_watch_dec_hpp
#define UZLMathLib_stop_watch_dec_hpp

UZLMATH_BEGIN

/*!
 * @brief   A measurment tool for getting execution time of code snippets.
 * @details The stop watch can be used to measure the execution time of a
 *          specified code snippet. The time can be returned in common time
 *          units like nano seconds, micro seconds, milli seconds, seconds,
 *          minutes and hours.
 *
 * @since   0.0.1
 *
 * @author  Denis-Michael Lux <denis.lux@icloud.com>
 * @date    01.05.15
 */
class stopwatch
{
    struct timeval start;   //!< The start time reference
    
    stopwatch();            
    
public:
    static stopwatch tic();
    ~stopwatch();
    
    double toc();
    double toc_micros();
    double toc_millis();
    double toc_seconds();
    double toc_minutes();
    double toc_hours();
};

UZLMATH_END

#endif
