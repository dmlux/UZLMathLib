//
//  fn_fourier_transforms.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 21.04.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_fn_fourier_transforms_hpp
#define uzlmath_fn_fourier_transforms_hpp

/*!
 * @brief       The fft namespace encapsulates a collection of fourier transforms for different
 *              datastructures like matrices and vectors.
 *
 * @since       0.0.1
 *
 * @author      Denis-Michael Lux <denis.lux@icloud.com>
 * @date        21.04.15
 */
namespace FourierTransforms
{
    /*!
     * @brief           Performing the FFT on a complex vector.
     * @details         Performing a FFT the given complex vector.
     *                  The FFT is defined as follows
     *                  \f{eqnarray*}{
     *                      Y_k = c\cdot \sum\limits_{j = 0}^{n-1}X_je^{-2\pi jk\sqrt{-1}/n}
     *                  \f}
     *                  where \f$Y\f$ is the Fourier transformed sequence of the input
     *                  input sequence \f$X\f$ of size \f$n\f$.
     *                  signal.
     *
     * @tparam          eT An element type which represents a number that provides all common
     *                  mathmatical operations.
     *
     * @param[in, out]  vec The complex vector that is supposed to be transformed.
     * @param[in]       scale The scaling factor for the matrix elements. Each element of
     *                  the matrix will be multiplied by the scaling factor after the FFT.
     *                  The value is 1 by default which means each value will be taken as
     *                  it is. Changing the factor to a custom factor \f$C\f$ means each
     *                  value will be multilpied by \f$C\f$.
     *
     * @sa              complex
     * @sa              vector
     *
     * @since           0.1.1
     *
     * @author          Denis-Michael Lux <denis.lux@icloud.com>
     * @date            23.06.15
     * 
     */
    template< typename eT >
    inline
    auto DFT(vector< eT >& vec, eT scale = eT(1, 0)) -> typename uzl_void_cx_num_only< eT >::result
    {
        if (vec.n_elements() == 0)
        {
            printf("** uzlmath error: vector for FFT has no size. **\n");
            return;
        }
        
        // create FFT data array
        double* inOut;
        
        // create array of given real data
        size_t i;
        if (is_double< eT >::value == true)
        {
            // If the POD type is double we can just cast the
            // complex array to an double because memory layout
            // is guaranteed by the compiler
            inOut = reinterpret_cast< double* >(vec.memptr());
        }
        else
        {
            // create array
            inOut = new double[2 * vec.n_elements()];
            
            // copy vector elements and store casted values
            for (i = 0; i < vec.n_elements(); ++i)
            {
                inOut[i * 2    ] = static_cast< double >(vec[i].re);
                inOut[i * 2 + 1] = static_cast< double >(vec[i].im);
            }
        }
        
        // call fftw function
        uzl_fftw_DFT(vec.n_elements(), inOut);
        
        // copy transformed values into given vector
        if (is_double< eT >::value == false)
        {
            for (i = 0; i < vec.n_elements(); ++i)
            {
                // fill vector with transformed values
                vec[i].re = inOut[i * 2    ];
                vec[i].im = inOut[i * 2 + 1];
            }
            
            // free allocated memory
            delete [] inOut;
        }
        
        if (scale.re != 1 || scale.im != 0)
        {
            // scale values
            for (i = 0; i < vec.n_elements(); ++i)
            {
                vec[i] *= scale;
            }
        }
    }
    
    /*!
     * @brief           Performing the FFT on a real vector.
     * @details         Performing an FFT the given real vector.
     *                  The IFFT is defined as follows
     *                  \f{eqnarray*}{
     *                      Y_k = c\cdot \sum\limits_{j = 0}^{n-1}X_je^{-2\pi jk\sqrt{-1}/n}
     *                  \f}
     *                  where \f$Y\f$ is the Fourier transformed sequence of the input
     *                  input sequence \f$X\f$ of size \f$n\f$.
     *                  signal.
     *
     * @tparam          eT An element type which represents a number that provides all common
     *                  mathmatical operations.
     *
     * @param[in, out]  vec The complex vector that is supposed to be transformed.
     * @param[in]       scale The scaling factor for the matrix elements. Each element of
     *                  the matrix will be multiplied by the scaling factor after the FFT.
     *                  The value is 1 by default which means each value will be taken as
     *                  it is. Changing the factor to a custom factor \f$C\f$ means each
     *                  value will be multilpied by \f$C\f$.
     *
     * @sa              complex
     * @sa              vector
     *
     * @since           0.1.1
     *
     * @author          Denis-Michael Lux <denis.lux@icloud.com>
     * @date            23.06.15
     *
     */
    template< typename eT >
    inline
    auto DFT(vector< eT >& vec, complex< double > scale = complex< double >(1, 0)) -> typename uzl_vec_cx_dbl_real_num_only< eT >::result
    {
        if (vec.n_elements() == 0)
        {
            printf("** uzlmath error: vector for FFT has no size. **\n");
            return vector< complex< double > >();
        }
        
        // create FFT data array
        vector< complex< double > > result(vec.n_elements(), vec.vecType());
        double* inOut = reinterpret_cast< double* >(result.memptr());
            
        // copy vector elements and store casted values
        size_t i;
        for (i = 0; i < vec.n_elements(); ++i)
        {
            inOut[i * 2    ] = static_cast< double >(vec[i]);
            inOut[i * 2 + 1] = static_cast< double >(0);
        }
        
        // call fftw function
        uzl_fftw_DFT(vec.n_elements(), inOut);
        
        // scale values
        if (scale.re != 1 || scale.im != 0)
        {
            // scale values
            for (i = 0; i < vec.n_elements(); ++i)
            {
                result[i] *= scale;
            }
        }
        
        // return result
        return result;
    }
    
    /*!
     * @brief           Performing the FFT2 on a complex matrix.
     * @details         Performing a FFT on each column and each row of the given 
     *                  matrix. The FFT on the row and column vectors is defined 
     *                  as follows
     *                  \f{eqnarray*}{
     *                      X(k) = \sum\limits_{j=1}^Nx(j)\cdot e^{-2\pi i(j-1)(k-1)/N}
     *                  \f}
     *                  where \f$X(k)\f$ is the Fourier transformed sequence and \f$x(j)\f$
     *                  the input sequence of size \f$N\f$
     *
     * @tparam          eT An element type which represents a number that provides all common
     *                  mathmatical operations.
     * @param[in, out]  mat The matrix that is supposed to be fourier transformed.
     * @param[in]       scale The scaling factor for the matrix elements. Each element of
     *                  the matrix will be multiplied by the scaling factor after the FFT.
     *                  The value is 1 by default which means each value will be taken as 
     *                  it is. Changing the factor to a custom factor \f$C\f$ means each 
     *                  value will be multilpied by \f$C\f$.
     *
     * @sa              complex
     * @sa              matrix
     *
     * @since           0.0.1
     *
     * @author          Denis-Michael Lux <denis.lux@icloud.com>
     * @date            21.04.15
     *
     * @todo            Direct cast complex matrix memory into double pointer
     */
    template< typename eT >
    inline
    auto DFT2(matrix< complex< eT > >& mat, complex< eT > scale = complex< eT >(1, 0)) -> typename uzl_void_real_num_only< eT >::result
    {
        // make fft array
        double fftInOut[mat.n_rows() * mat.n_cols() * 2];
        
        // extract complex values
        size_t i;
        for (i = 0; i < mat.n_rows() * mat.n_cols(); ++i)
        {
            fftInOut[i * 2]     = mat.memptr()[i].re;
            fftInOut[i * 2 + 1] = mat.memptr()[i].im;
        }
        
        uzl_fftw_DFT2(mat.n_cols(), mat.n_rows(), fftInOut);
        
        // translate back into complex matrix
        for (i = 0; i < mat.n_rows() * mat.n_cols(); ++i)
        {
            complex< eT > comp(fftInOut[i * 2], fftInOut[i * 2 + 1]);
            mat.memptr()[i] = scale * comp;
        }
    }
    
    /*!
     * @brief           Performing the IFFT2 on a complex matrix.
     * @details         Performing a IFFT on each column and each row of the given
     *                  matrix. The FFT on the row and column vectors is defined
     *                  as follows
     *                  \f{eqnarray*}{
     *                      x(k) = \frac{1}{N}\sum\limits_{j=1}^NX(j)\cdot e^{2\pi i(j-1)(k-1)/N}
     *                  \f}
     *                  where \f$x(k)\f$ is the input sequence and \f$X(j)\f$
     *                  the Fourier transformed sequence of size \f$N\f$
     *
     * @tparam          eT An element type which represents a number that provides all common
     *                  mathmatical operations.
     * @param[in, out]  mat The matrix that is supposed to be inverse fourier transformed.
     * @param[in]       scale The scaling factor for the matrix elements. Each element of
     *                  the matrix will be multiplied by the scaling factor after the IFFT.
     *                  The value is 1 by default which means each value will be multiplied
     *                  by \f$\frac{1}{N}\f$. Changing the factor to a custom factor \f$C\f$
     *                  means each value will be multilpied by \f$\frac{C}{N}\f$.
     *
     * @sa              complex
     * @sa              matrix<>
     *
     * @since           0.0.1
     *
     * @author          Denis-Michael Lux <denis.lux@icloud.com>
     * @date            21.04.15
     *
     * @todo            Direct cast complex matrix memory into double pointer
     */
    template< typename eT >
    inline
    auto IDFT2(matrix< complex< eT > >& mat, complex< eT > scale = complex< eT >(1,0)) -> typename uzl_void_real_num_only< eT >::result
    {
        // make fft array
        double fftInOut[mat.n_rows() * mat.n_cols() * 2];
        
        // extract complex values
        size_t i;
        for (i = 0; i < mat.n_rows() * mat.n_cols(); ++i)
        {
            fftInOut[i * 2]     = mat.memptr()[i].re;
            fftInOut[i * 2 + 1] = mat.memptr()[i].im;
        }
        
        uzl_fftw_IDFT2(mat.n_cols(), mat.n_rows(), fftInOut);
        
        // translate back into complex matrix
        for (i = 0; i < mat.n_rows() * mat.n_cols(); ++i)
        {
            complex< eT > comp(fftInOut[i * 2], fftInOut[i * 2 + 1]);
            mat.memptr()[i] = scale * complex< eT >(1.0 / (mat.n_rows() * mat.n_cols()), 0) * comp;
        }
    }
    
    // Forward fast Fourier transform on SO(3)
    auto SOFT(grid3D< complex< double > > sample, SOFTFourierCoefficients& fc) -> void;
    
    // Inverse fast Fourier transform on SO(3)
    auto ISOFT(const SOFTFourierCoefficients& fc, grid3D< complex< double > >& synthesis) -> void;
}

#endif
