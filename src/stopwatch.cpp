//
//  stopwatch.cpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 18.06.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_stopwatch_cpp
#define UZLMathLib_stopwatch_cpp

#include <uzlmath>

UZLMATH_BEGIN
    
/*!
 * @brief           Private constructor for the stopwatch class
 * @details         Constructs the stopwatch by setting the start
 *                  time to the current time.
 */
stopwatch::stopwatch()
{
    gettimeofday(&start, NULL);
}

/*!
 * @brief           The stopwatch::tic() method creates a stopwatch object
 * @details         Constructs an stopwatch object ans sets the start time
 *                  to current timestamp.
 *
 * @return          A stopwatch object.
 */
stopwatch stopwatch::tic()
{
    // create object
    return stopwatch();
}

/*!
 * @brief           The default destructor method.
 */
stopwatch::~stopwatch()
{}

/*!
 * @brief           Stops the time in seconds from the stopwatch::tic() or the last
 *                  toc_<X>() time
 * @details         Stops the time and calculates the difference from toc or last tic 
 *                  command in seconds.
 *
 * @return          The seconds elapsed from stopwatch::tic() or other toc_<X>() functions
 */
double stopwatch::toc()
{
    // set end time
    struct timeval end;
    gettimeofday(&end, NULL);
    
    // calculate difference
    long elapsed = ((end.tv_sec - start.tv_sec) * 1000000L + end.tv_usec) - start.tv_usec;
    
    // reset time
    start = end;
    
    // return result
    return elapsed / 1e+6;
}

/*!
 * @brief           Stops the time in microseconds from the stopwatch::tic() or the 
 *                  last toc_<X>() time
 * @details         Stops the time and calculates the difference from toc or last tic
 *                  command in microseconds.
 *
 * @return          The microseconds elapsed from stopwatch::tic() or other toc_<X>() functions
 */
double stopwatch::toc_micros()
{
    // set end time
    struct timeval end;
    gettimeofday(&end, NULL);
    
    // calculate difference
    long elapsed = ((end.tv_sec - start.tv_sec) * 1000000L + end.tv_usec) - start.tv_usec;
    
    // reset time
    start = end;
    
    // return result
    return elapsed;
}

/*!
 * @brief           Stops the time in milliseconds from the stopwatch::tic() or the
 *                  last toc_<X>() time
 * @details         Stops the time and calculates the difference from toc or last tic
 *                  command in milliseconds.
 *
 * @return          The milliseconds elapsed from stopwatch::tic() or other toc_<X>()
 */
double stopwatch::toc_millis()
{
    // set end time
    struct timeval end;
    gettimeofday(&end, NULL);
    
    // calculate difference
    long elapsed = ((end.tv_sec - start.tv_sec) * 1000000L + end.tv_usec) - start.tv_usec;
    
    // reset time
    start = end;
    
    // return result
    return elapsed / 1e+3;
}

/*!
 * @brief           Stops the time in seconds from the stopwatch::tic() or the last
 *                  toc_<X>() time
 * @details         Stops the time and calculates the difference from toc or last tic
 *                  command in seconds.
 *
 * @return          The seconds elapsed from stopwatch::tic() or other toc_<X>() functions
 */
double stopwatch::toc_seconds()
{
    // set end time
    struct timeval end;
    gettimeofday(&end, NULL);
    
    // calculate difference
    long elapsed = ((end.tv_sec - start.tv_sec) * 1000000L + end.tv_usec) - start.tv_usec;
    
    // reset time
    start = end;
    
    // return result
    return elapsed / 1e+6;
}

/*!
 * @brief           Stops the time in minutes from the stopwatch::tic() or the last
 *                  toc_<X>() time
 * @details         Stops the time and calculates the difference from toc or last tic
 *                  command in minutes.
 *
 * @return          The minutes elapsed from stopwatch::tic() or other toc_<X>() functions
 */
double stopwatch::toc_minutes()
{
    // set end time
    struct timeval end;
    gettimeofday(&end, NULL);
    
    // calculate difference
    long elapsed = ((end.tv_sec - start.tv_sec) * 1000000L + end.tv_usec) - start.tv_usec;
    
    // reset time
    start = end;
    
    // return result
    return elapsed / 6e+7;
}

/*!
 * @brief           Stops the time in hours from the stopwatch::tic() or the last
 *                  toc_<X>() time
 * @details         Stops the time and calculates the difference from toc or last tic
 *                  command in minutes.
 *
 * @return          The minutes elapsed from stopwatch::tic() or other toc_<X>() functions
 */
double stopwatch::toc_hours()
{
    // set end time
    struct timeval end;
    gettimeofday(&end, NULL);
    
    // calculate difference
    long elapsed = ((end.tv_sec - start.tv_sec) * 1000000L + end.tv_usec) - start.tv_usec;
    
    // reset time
    start = end;
    
    // return result
    return elapsed / 36e+8;
}
    
UZLMATH_END

#endif
