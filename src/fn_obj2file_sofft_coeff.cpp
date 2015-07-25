//
//  fn_obj2file_soft_coeff.cpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 18.06.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_fn_obj2file_soft_coeff_cpp
#define UZLMathLib_fn_obj2file_soft_coeff_cpp

#include <uzlmath>

UZLMATH_BEGIN
    
/*!
 * @brief           Writes the Fourier coefficients container to disk.
 * @details         Writes the contents of the Fourier coefficients container
 *                  to a file with the given name. 
 *
 * @param[in]       coef The Fourier coefficients container containing the Fourier
 *                  coefficients of the SOFT for the given bandwidth.
 * @param[in]       fileName The file name where the coefficients are getting written
 *                  to.
 *
 * @ingroup         SOFTFourierCoefficients
 */
auto obj2file(const SOFTFourierCoefficients& coef, const std::string& fileName) -> void
{
    // create file and open it with write flag
    FILE* fp = fopen(fileName.c_str(), "w");
    
    // print labels into  file labels
    fprintf(fp, "-----------------------------------------------------------------\n");
    fprintf(fp, "  l    M    M'     Fourier coefficient\n");
    fprintf(fp, "-----------------------------------------------------------------\n");
    
    // iterate over grid starting with layers
    int l, M, Mp;
    for (l = 0; l < coef.bandwidth; ++l)
    {
        for (M = -l; M <= l; ++M)
        {
            for (Mp = -l; Mp <= l; ++Mp)
            {
                // print orders and degree of coefficient
                fprintf(fp, "%4d %4d %4d     ", l, M, Mp);
                // print space if real part of coefficient is positive
                fprintf(fp, "%s", (coef(l, M, Mp).re > 0 ? " " : ""));
                // print coefficient
                fprintf(fp, "%.16f%s%.16fi\n", coef(l, M, Mp).re, (coef(l, M, Mp).im >= 0 ? " + " : " - "), fabs(coef(l, M, Mp).im));
            }
        }
    }
    
    // close file and everything is done
    fclose(fp);
}

UZLMATH_END
    
#endif
