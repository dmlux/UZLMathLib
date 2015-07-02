//
//  soft_fourier_coefficients.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 05.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_SOFTFourierCoefficients_hpp
#define uzlmath_SOFTFourierCoefficients_hpp

/*!
 * @brief       Collection of functions and classes for SOFT Fourier coefficients
 *              container
 * @defgroup    SOFTFourierCoefficients SOFT Fourier coefficients container
 * @{
 */

/*!
 * @brief   A datastructure to manage fourier coefficients that 
 *          where produced by the SOFT algorithm described by 
 *          P. J. Kostelec and D. N. Rockmore in the paper
 *          'FFTs on the Rotation Group'
 * @details This class provides a manager class that organizes
 *          and stores the fourier coefficients that where produced
 *          by the SOFT algorithm. The fourier coefficients are
 *          indexed over three parameter.
 *
 * @since   0.0.1
 *
 * @author  Denis-Michael Lux <denis.lux@icloud.com>
 * @date    05.05.15
 */
struct SOFTFourierCoefficients
{
private:
    const int max_l;                  //!< Maximum l index
    const int B;                      //!< Bandwidth of function
    matrix< complex< double > >* mem; //!< Index space
    
public:
                             SOFTFourierCoefficients();
                             SOFTFourierCoefficients(int bandlimit);
                            ~SOFTFourierCoefficients();
    
          complex< double >& operator()(const int& l, const int& M, const int& Mp);
    const complex< double >& operator()(const int& l, const int& M, const int& Mp) const;
    
    const int&               bandwidth() const;
    
    friend std::ostream&     operator<<(std::ostream& o, const SOFTFourierCoefficients& fc);
};

/*!
 * @}
 */

#endif
