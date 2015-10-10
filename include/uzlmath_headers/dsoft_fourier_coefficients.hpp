//
//  soft_fourier_coefficients.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 05.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_dsoft_fourier_coefficients_hpp
#define UZLMathLib_dsoft_fourier_coefficients_hpp

UZLMATH_BEGIN

/*!
 * @brief       Collection of functions and classes for SOFT Fourier coefficients
 *              container
 * @defgroup    DSOFTFourierCoefficients SOFT Fourier coefficients container
 * @{
 */

/*!
 * @brief       A datastructure to manage fourier coefficients that
 *              where produced by the SOFT algorithm described by
 *              P. J. Kostelec and D. N. Rockmore in the paper
 *              'FFTs on the Rotation Group'
 * @details     This class provides a manager class that organizes
 *              and stores the fourier coefficients that where produced
 *              by the SOFT algorithm. The fourier coefficients are
 *              indexed over three parameter.
 *
 * @since       0.0.1
 *
 * @author      Denis-Michael Lux <denis.lux@icloud.com>
 * @date        05.05.15
 */
struct DSOFTFourierCoefficients
{
private:
    matrix< complex< double > >* mem; //!< Coefficients storage
    
public:
    // public ivars
    const int bandwidth;              //!< Bandwidth of function
    
    // constructors
    DSOFTFourierCoefficients();
    DSOFTFourierCoefficients(int bandlimit);
    
    // destructor
    ~DSOFTFourierCoefficients();
    
    // methods
          complex< double >&            operator()(const int& l, const int& M, const int& Mp);
    const complex< double >&            operator()(const int& l, const int& M, const int& Mp) const;
    
    friend std::ostream& operator<<(std::ostream& o, const DSOFTFourierCoefficients& fc);
};

// prototype for the overloaded stream operator
std::ostream& operator<<(std::ostream& o, const DSOFTFourierCoefficients& fc);

/*!
 * @}
 */

UZLMATH_END

#endif
