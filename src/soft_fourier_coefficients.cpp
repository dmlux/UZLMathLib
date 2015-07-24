//
//  soft_fourier_coefficients.cpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 18.06.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_SOFTFourierCoefficients_cpp
#define uzlmath_SOFTFourierCoefficients_cpp

#include <uzlmath>

UZLMATH_BEGIN

/*!
 * @brief           Default constructor
 * @details         Constructs a Fourier coefficients container that is
 *                  empty
 */
SOFTFourierCoefficients::SOFTFourierCoefficients()
    : mem(nullptr)
    , bandwidth(0)
{}


/*!
 * @brief           Constructor for a SOFTFourierCoefficients container
 * @details         Creating and allocating memory for Fourier coefficents
 *                  for a SOFT which is described in the paper "FFT's on
 *                  the rotation group". This container is memory manager
 *                  for \f$\hat{f}^l_{M,M'}\f$.
 *
 * @param[in]       bandlimit The bandlimit of the function which coefficients
 *                  are supposed to be stored in this coefficient container.
 */
SOFTFourierCoefficients::SOFTFourierCoefficients(int bandlimit)
    : bandwidth(bandlimit)
{
    mem  = new matrix< complex< double > >[bandlimit];
    
    for (int i = 0; i < bandlimit; ++i)
    {
        mem[i] = matrix< complex< double > >(2 * i + 1, 2 * i + 1);
    }
}

/*!
 * @brief           Destructor for the SOFTFourierCoefficients manager
 * @details         Frees the memory that is allocated for the coefficents.
 */
SOFTFourierCoefficients::~SOFTFourierCoefficients()
{
    delete [] mem;
}

/*!
 * @brief           Accessor operator for the SOFTFourierCoefficients manager
 * @details         Makes the memory for the coefficents accessable by using
 *                  degree and orders.
 *
 * @param[in]       l Degree of the fourier coefficent
 * @param[in]       M First oder of the fourier coefficient
 * @param[in]       Mp Second order of the fourier coefficient
 *
 * @return          Pointer to memory of fourier coefficient with degree \f$l\f$
 *                  and orders \f$M\f$ and \f$M'\f$
 */
complex< double >& SOFTFourierCoefficients::operator()(const int& l, const int& M, const int& Mp)
{
    if (M > l || Mp > l || M < -l || Mp < -l)
    {
        uzlmath_error("%s", "illegal parameter configuration for SOFTFourierCoefficients access. M > l, M < -l, Mp < -l or Mp > l.");
    }
    
    // if M or Mp are negative count from behind
    size_t idx_M  = (M  >= 0 ? M  : mem[l].n_rows() + M );
    size_t idx_Mp = (Mp >= 0 ? Mp : mem[l].n_cols() + Mp);
    
    return mem[l](idx_M, idx_Mp);
}

/*!
 * @brief           Accessor operator for the SOFTFourierCoefficients manager
 * @details         Makes the memory for the coefficents readable by using
 *                  degree and orders.
 *
 * @param[in]       l Degree of the fourier coefficent
 * @param[in]       M First oder of the fourier coefficient
 * @param[in]       Mp Second order of the fourier coefficient
 *
 * @return          The value of fourier coefficient with degree \f$l\f$
 *                  and orders \f$M\f$ and \f$M'\f$
 */
const complex< double >& SOFTFourierCoefficients::operator()(const int& l, const int& M, const int& Mp) const
{
    if (M > l || Mp > l || M < -l || Mp < -l)
    {
        uzlmath_error("%s", "illegal parameter configuration for SOFTFourierCoefficients access. M > l, M < -l, Mp < -l or Mp > l.");
    }
    
    // if M or Mp are negative count from behind
    size_t idx_M  = (M  >= 0 ? M  : mem[l].n_rows() + M );
    size_t idx_Mp = (Mp >= 0 ? Mp : mem[l].n_cols() + Mp);
    
    return mem[l](idx_M, idx_Mp);
}

/*!
 * @brief           Outstream operator overload for SOFTFourierCoefficients.
 * @details         The out-steam operator is used to print the coefficents in
 *                  a nice form over the std::cout stream.
 * 
 * @param[in,out]   o The stream object
 * @param[in]       fc The SOFTFourierCoefficients manager
 *
 * @return          The reference to the given out-stream.
 */
std::ostream& operator<<(std::ostream& o, const SOFTFourierCoefficients& fc)
{
    std::ios::fmtflags f( std::cout.flags() );
    o << std::endl << std::setprecision(4);
    
    for (int i = 0; i < fc.bandwidth; ++i)
    {
        o << "SOFTFourierCoefficients[M_{0,1,2,...,-2,-1} x M'_{0,1,2,...,-2,-1}] ~> [l = " << i << "]" << std::endl;
        o << fc.mem[i] << std::endl;
    }
    
    std::cout.flags( f );
    return o;
}

UZLMATH_END

#endif