//
//  fn_fourier_transforms.cpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 18.06.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_fn_fourier_transforms_cpp
#define UZLMathLib_fn_fourier_transforms_cpp

#include <uzlmath>

UZLMATH_NAMESPACE(FourierTransforms)
        
/*!
 * @brief           The SOFT (<b>S0</b>(3) <b>F</b>ourier <b>T</b>ransform)
 *                  describes the FFT on the rotation group \f$\mathcal{SO}(3)\f$
 * @details         The method to compute the SOFT on the rotation group is described in
 *                  detail in the paper 'FFTs on the Rotation Group' written by Peter J.
 *                  Kostelec and Daniel N. Rockmore. The implementation that is underneath
 *                  this function is
 *                  \f[
 *                      \hat{f}^l_{M,M'} = \frac{\pi}{(2B)^2}\sum\limits_{k = 0}^{2B-1}w_B(k)
 *                          \tilde{d}^l_{M,M'}(\beta_k)\sum\limits_{j2 = 0}^{2B-1}e^{iM'\gamma_{j_2}}
 *                          \sum\limits_{j_1 = 0}^{2B-1}e^{iM\alpha_{j_1}}f(\alpha_{j_1},\beta_k,\gamma_{j_2})
 *                  \f]
 *                  where \f$B\f$ is the bandwidth of function \f$f(\alpha_{j_1},\beta_k,\gamma_{j_2})\f$.
 *                  The number of cofficients \f$\hat{f}^l_{M,M'}\f$ can be calculated by
 *                  \f[
 *                      |\{\hat{f}^l_{M,M'}\}_{l\in\{0,\dots,B-1\}, M,M'\in\{-l,\dots,l\}}| = \sum\limits_{i = 0}^{B-1}(2 * i + 1)^2
 *                  \f]
 *                  The implementation itself uses the symmetry properties of the wigner
 *                  d-function to reduce the number of wigner d-function evaluations. The
 *                  symmetries that are used are
 *                  \f{eqnarray*}{
 *                      d^{J}_{MM'}(\beta) &=& (-1)^{M-M'}d^J_{-M-M'}(\beta)\\
 *                      &=& (-1)^{M-M'}d^J_{M'M}(\beta)\\
 *                      &=& d^J_{-M'-M}(\beta)\\
 *                      &=& (-1)^{J-M'}d^J_{-MM'}(\pi-\beta)\\
 *                      &=& (-1)^{J+M}d^J_{M-M'}(\pi-\beta)\\
 *                      &=& (-1)^{J-M'}d^J_{-M'M}(\pi-\beta)\\
 *                      &=& (-1)^{J+M}d^J_{M'-M}(\pi-\beta)
 *                  \f}
 *
 * @param[in]       sample A discrete sample of function \f$f\f$ which has the
 *                  dimension of \f$2B\times 2B\times 2B\f$.
 * @param[out]      fc A Fourier coefficent managment container with capacaty for
 *                  all Fourier coefficients of \f$f\f$.
 *
 * @sa              DWT::quadrature_weights
 * @sa              DWT::wigner_d_matrix
 * @sa              SOFTFourierCoefficients
 * @sa              complex
 * @sa              matrix
 * @sa              grid3D
 *
 * @ingroup         fourierTransforms
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            14.05.2015
 *
 * @since           0.0.1
 */
auto SOFT(grid3D< complex< double > > sample, SOFTFourierCoefficients& fc, int threads) -> void
{
    /*****************************************************************
     ** Check parameters                                            **
     *****************************************************************/
    // Check if the grid has same size in each dimension
    if (sample.rows != sample.cols || sample.rows != sample.lays)
    {
        uzlmath_warning("%s", "all SOFT sample grid dimensions should be equal.");
        return;
    }
    
    // Check if grid has odd dimensions
    if (sample.rows & 1)
    {
        uzlmath_warning("%s", "SOFT sample grid dimensions are not even.");
        return;
    }
    
    // Extract bandwidth
    int bandwidth = static_cast< int >(sample.cols / 2);
    
    // Check if Fourier coefficients container dimension matches sample dimension
    if (bandwidth != fc.bandwidth)
    {
        uzlmath_warning("%s", "SOFT Fourier coefficients container bandwidth does not match to sample grid bandwidth.");
        return;
    }
    
    // print warinings for serial implementation
    #ifndef _OPENMP
    if (threads != 1)
    {
        uzlmath_warning("%s", "compiler does not support OpenMP. Changing the number of threads for the SOFT has no effect.");
    }
    #endif
    
    /*****************************************************************
     ** FFT2 transform layers of sample grid for fixed k            **
     *****************************************************************/
    sample.layer_wise_DFT2();
    
    /*****************************************************************
     ** M = 0, M' = 0                                               **
     *****************************************************************/
    vector< double > weights = DWT::quadrature_weights(bandwidth);
    matrix< double >      dw = DWT::weighted_wigner_d_matrix(bandwidth, 0, 0, weights) * -1;
    vector< complex< double > > s(2 * bandwidth, vec_type::COLUMN);
    
    // defining norm factor
    complex< double > norm(M_PI / (2 * bandwidth * bandwidth), 0);
    
    // defining needed indices
    int MMp, i, M, Mp;
    
    // precompute the double bandwidth
    const int bw2 = 2 * bandwidth;
    
    // DWT for M = 0, M' = 0
    for (const complex< double >* e = s.mem; e != s.mem + bw2; ++e) { access::rw(*e) = sample(0, 0, e - s.mem);             }
    vector< complex< double > > sh = dw * s;
    for (i = 1; i <= sh.size; ++i)                                  { fc(bandwidth - i, 0, 0) = norm * sh[sh.size - i];     }
    
    /*****************************************************************
     ** Iterate over all combinations of M and M'                   **
     *****************************************************************/
    #pragma omp parallel default(shared) if(bandwidth >= SOFT_THRESHOLD) num_threads(threads)
    {
        
        #pragma omp for private(i, M, dw, sh) firstprivate(s) schedule(dynamic) nowait
        for (M = 1; M < bandwidth; ++M)
        {
            dw = DWT::weighted_wigner_d_matrix(bandwidth, M, 0, weights) * -1;
            
            /*****************************************************************
             ** Make use of symmetries                                      **
             *****************************************************************/
            // case f_{M,0}
            for (const complex< double >* e = s.mem; e != s.mem + bw2; ++e) { access::rw(*e) = sample(0, M, e - s.mem);                 }
            sh = dw * s;
            for (i = 1; i <= sh.size; ++i)                                  { fc(bandwidth - i, M, 0) = norm * sh[sh.size - i];         }
            
            // case f_{0,M}
            for (const complex< double >* e = s.mem; e != s.mem + bw2; ++e) { access::rw(*e) = sample(M, 0, e - s.mem);                 }
            sh = dw * s;
            if  (M & 1)                                                     { sh *= -1;                                                 }
            for (i = 1; i <= sh.size; ++i)                                  { fc(bandwidth - i, 0, M) = norm * sh[sh.size - i];         }
            
            // case f_{-M,0}
            fliplr(dw);
            for (const complex< double >* e = s.mem; e != s.mem + bw2; ++e) { access::rw(*e) = sample(0, bw2 - M, e - s.mem);           }
            sh = dw * s;
            if (M & 1)
            {
                for (i = 0; i < sh.size; i += 2)                            { sh[i] *= -1;                                              }
            }
            else
            {
                for (i = 1; i < sh.size; i += 2)                            { sh[i] *= -1;                                              }
            }
            for (i = 1; i <= sh.size; ++i)                                  { fc(bandwidth - i, -M, 0) = norm * sh[sh.size - i];        }
            
            // case f_{0,-M}
            for (const complex< double >* e = s.mem; e != s.mem + bw2; ++e) { access::rw(*e) = sample(bw2 - M, 0, e - s.mem);           }
            sh = dw * s;
            for (i = 1; i < sh.size; i +=2)                                 { sh[i] *= -1;                                              }
            for (i = 1; i <= sh.size; ++i)                                  { fc(bandwidth - i, 0, -M) = norm * sh[sh.size - i];        }
            
            // get new wigner matrix
            dw = DWT::weighted_wigner_d_matrix(bandwidth, M, M, weights) * -1;
            
            // case f_{M, M}
            for (const complex< double >* e = s.mem; e != s.mem + bw2; ++e) { access::rw(*e) = sample(M, M, e - s.mem);                 }
            sh = dw * s;
            for (i = 1; i <= sh.size; ++i)                                  { fc(bandwidth - i, M, M) = norm * sh[sh.size - i];         }
            
            // case f_{-M, -M}
            for (const complex< double >* e = s.mem; e != s.mem + bw2; ++e) { access::rw(*e) = sample(bw2 - M, bw2 - M, e - s.mem);     }
            sh = dw * s;
            for (i = 1; i <= sh.size; ++i)                                  { fc(bandwidth - i, -M, -M) = norm * sh[sh.size - i];       }
            
            // Modify dw for the last two cases. flip matrix from left to right and negates sign of
            // every second row with odd row indices.
            fliplr_ne2ndorow(dw);
            
            // A little arithmetic error is occuring in the following calculation... I do not exactly know why
            // case f_{M, -M}
            for (const complex< double >* e = s.mem; e != s.mem + bw2; ++e) { access::rw(*e) = sample(bw2 - M, M, e - s.mem);           }
            sh = dw * s;
            for (i = 1; i <= sh.size; ++i)                                  { fc(bandwidth - i, M, -M) = norm * sh[sh.size - i];        }
            
            // case f_{-M, M}
            for (const complex< double >* e = s.mem; e != s.mem + bw2; ++e) { access::rw(*e) = sample(M, bw2 - M, e - s.mem);           }
            sh = dw * s;
            for (i = 1; i <= sh.size; ++i)                                  { fc(bandwidth - i, -M, M) = norm * sh[sh.size - i];        }
        }
        
        // Fused two loops per hand
        //
        // for (M = 1; M < bandwidth; ++M)
        //     for (Mp = 1; Mp < M; ++Mp)
        //
        // which now is equivalent to the following loop
        #pragma omp for private(i, MMp, M, Mp, dw, sh) schedule(dynamic) firstprivate(s) nowait
        for (MMp = 0; MMp < (bandwidth - 2) * (bandwidth - 1) / 2; ++MMp)
        {
            // reconstructing nested loop indices
            int i = MMp / (bandwidth - 1) + 1;
            int j = MMp % (bandwidth - 1) + 1;
            
            // get M and M'
            M  = j > i ? bandwidth - i : i + 1;
            Mp = j > i ? bandwidth - j : j    ;
            
            // get new wigner d-matrix
            dw = DWT::weighted_wigner_d_matrix(bandwidth, M, Mp, weights);
            
            // case f_{M, Mp}
            for (const complex< double >* e = s.mem; e != s.mem + bw2; ++e) { access::rw(*e) = sample(Mp, M, e - s.mem);                        }
            sh  = dw * s;
            sh *= -1;
            for (i = 1; i <= sh.size; ++i)                                  { fc(bandwidth - i, M, Mp) = norm * sh[sh.size - i];                }
            
            // case f_{Mp, M}
            for (const complex< double >* e = s.mem; e != s.mem + bw2; ++e) { access::rw(*e) = sample(M, Mp, e - s.mem);                        }
            sh = dw * s;
            if  (!((M - Mp) & 1))                                           { sh *= -1;                                                         }
            for (i = 1; i <= sh.size; ++i)                                  { fc(bandwidth - i, Mp, M) = norm * sh[sh.size -  i];               }
            
            // case f_{-M, -Mp}
            for (const complex< double >* e = s.mem; e != s.mem + bw2; ++e) { access::rw(*e) = sample(bw2 - Mp, bw2 - M, e - s.mem);            }
            sh = dw * s;
            if  (!((M - Mp) & 1))                                           { sh *= -1;                                                         }
            for (i = 1; i <= sh.size; ++i)                                  { fc(bandwidth - i, -M, -Mp) = norm * sh[sh.size - i];              }
            
            // case f_{-Mp, -M}
            for (const complex< double >* e = s.mem; e != s.mem + bw2; ++e) { access::rw(*e) = sample(bw2 - M, bw2 - Mp, e - s.mem);            }
            sh  = dw * s;
            sh *= -1;
            for (i = 1; i <= sh.size; ++i)                                  { fc(bandwidth - i, -Mp, -M) = norm * sh[sh.size - i];              }
            
            // modify wigner d-matrix for next four cases. This just works because the weight
            // function is also symmetric like the wigner-d matrix. flip left-right the dw
            // matrix and negate each even value with even row index.
            fliplr_ne2nderow(dw);
            
            // case f_{Mp, -M}
            for (const complex< double >* e = s.mem; e != s.mem + bw2; ++e) { access::rw(*e) = sample(bw2 - M, Mp, e - s.mem);                  }
            sh = dw * s;
            for (i = 1; i <= sh.size; ++i)                                  { fc(bandwidth - i, Mp, -M) = norm * sh[sh.size - i];               }
            
            // case f_{M, -Mp}
            for (const complex< double >* e = s.mem; e != s.mem + bw2; ++e) { access::rw(*e) = sample(bw2 - Mp, M, e - s.mem);                  }
            sh = dw * s;
            for (i = 1; i <= sh.size; ++i)                                  { fc(bandwidth - i, M, -Mp) = norm * sh[sh.size - i];               }
            
            // alter signs
            if ((M - Mp) & 1)
            {
                for (i = 0; i < dw.rows * dw.cols; ++i)
                {
                    access::rw(dw.mem[i]) *= -1;
                }
            }
            
            // case f_{-Mp, M}
            for (const complex< double >* e = s.mem; e != s.mem + bw2; ++e) { access::rw(*e) = sample(M, bw2 - Mp, e - s.mem);                  }
            sh = dw * s;
            for (i = 1; i <= sh.size; ++i)                                  { fc(bandwidth - i, -Mp, M) = norm * sh[sh.size - i];               }
            
            // case f_{-M, Mp}
            for (const complex< double >* e = s.mem; e != s.mem + bw2; ++e) { access::rw(*e) = sample(Mp, bw2 - M, e - s.mem);                  }
            sh = dw * s;
            for (i = 1; i <= sh.size; ++i)                                  { fc(bandwidth - i, -M, Mp) = norm * sh[sh.size - i];               }
        }
    }
}

/*!
 * @brief           The inverse SOFT (<b>S0</b>(3) <b>F</b>ourier <b>T</b>ransform)
 *                  describes the inverse FFT on the rotation group \f$\mathcal{SO}(3)\f$
 * @details         The method to compute the inverse SOFT on the rotation group is described
 *                  in detail in the paper 'FFTs on the Rotation Group' written by Peter J.
 *                  Kostelec and Daniel N. Rockmore. The implementation that is underneath
 *                  this function is
 *                  \f[
 *                      f(\alpha,\beta,\gamma) = \sum\limits_{J\geq 0}\sum\limits^J_{M=-J}\sum\limits^J_{M'=-J}
 *                          \hat{f}^J_{MM'}\tilde{D}^J_{MM'}(\alpha,\beta,\gamma)
 *                  \f]
 *                  where \f$\hat{f}^J_{MM'}\f$ is the Fourier coefficient of degree
 *                  \f$J\f$ and orders \f$M,\;M'\f$. For more detailed information about the
 *                  Fourier coefficients read the documentation of spherical_harmonics::SOFT.
 *                  The implementation itself uses the symmetry properties of the wigner
 *                  d-function to reduce the number of wigner d-function evaluations. The
 *                  symmetries that are used are
 *                  \f{eqnarray*}{
 *                      d^{J}_{MM'}(\beta) &=& (-1)^{M-M'}d^J_{-M-M'}(\beta)\\
 *                      &=& (-1)^{M-M'}d^J_{M'M}(\beta)\\
 *                      &=& d^J_{-M'-M}(\beta)\\
 *                      &=& (-1)^{J-M'}d^J_{-MM'}(\pi-\beta)\\
 *                      &=& (-1)^{J+M}d^J_{M-M'}(\pi-\beta)\\
 *                      &=& (-1)^{J-M'}d^J_{-M'M}(\pi-\beta)\\
 *                      &=& (-1)^{J+M}d^J_{M'-M}(\pi-\beta)
 *                  \f}
 *
 * @param[in]       fc A Fourier coefficent managment container with all Fourier coefficients
 *                  of the SOFT.
 * @param[out]      synthesis The synthesized sample for the given Fourier coefficients.
 *
 * @sa              DWT::wigner_d_matrix
 * @sa              SOFTFourierCoefficients
 * @sa              FourierTransforms::SOFT
 * @sa              complex
 * @sa              matrix
 * @sa              grid3D
 *
 * @since           0.0.1
 *
 * @ingroup         fourierTransforms
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            23.05.2015
 */
auto ISOFT(const SOFTFourierCoefficients& fc, grid3D< complex< double > >& synthesis, int threads) -> void
{
    /*****************************************************************
     ** Check parameters                                            **
     *****************************************************************/
    // Check if the grid has same size in each dimension
    if (synthesis.rows != synthesis.cols || synthesis.rows != synthesis.lays)
    {
        uzlmath_warning("%s", "all ISOFT synthesis grid dimensions should be equal.");
        return;
    }
    
    // Check if grid has odd dimensions
    if (synthesis.rows & 1)
    {
        uzlmath_warning("%s", "ISOFT synthesis grid dimensions are not even.");
        return;
    }
    
    // Extract bandwidth
    int bandwidth = static_cast< int >(synthesis.cols / 2);
    
    // Check if Fourier coefficients container dimension matches sample dimension
    if (bandwidth != fc.bandwidth)
    {
        uzlmath_warning("%s", "ISOFT Fourier coefficients container bandwidth does not match to synthesis grid bandwidth.");
        return;
    }
    
    // print warinings for serial implementation
    #ifndef _OPENMP
    if (threads != 1)
    {
        uzlmath_warning("%s", "compiler does not support OpenMP. Changing the number of threads for the ISOFT has no effect.");
    }
    #endif
    
    /*****************************************************************
     ** M = 0, M' = 0                                               **
     *****************************************************************/
    matrix< double > d = DWT::wigner_d_matrix(bandwidth, 0, 0) * -1;
    vector< complex< double > > sh(d.rows, vec_type::COLUMN);
    
    d.transpose();
    
    // defining norm factor
    complex< double > norm((2 * bandwidth * bandwidth) / M_PI, 0);
    
    // defining needed indices
    int MMp, i, M, Mp;
    
    // inverse DWT for M = 0, M' = 0
    for (i = 1; i <= sh.size; ++i)      { sh[sh.size - i] = norm * fc(bandwidth - i, 0, 0);                     }
    vector< complex< double > > s = d * sh;
    for (i = 0; i < 2 * bandwidth; ++i) { synthesis(0, 0, i) = s[i];                                            }
    
    /*****************************************************************
     ** Iterate over all combinations of M and M'                   **
     *****************************************************************/
    #pragma omp parallel default(shared) if(bandwidth >= SOFT_THRESHOLD) num_threads(threads)
    {
        
        #pragma omp for private(i, M, d, s, sh) schedule(dynamic) nowait
        for (M = 1; M < bandwidth; ++M)
        {
            d  = DWT::wigner_d_matrix(bandwidth, M, 0) * -1;
            sh = vector< complex< double > >(d.rows, vec_type::COLUMN);
            d.transpose();
            
            /*****************************************************************
             ** Make use of symmetries                                      **
             *****************************************************************/
            // case f_{M,0}
            for (i = 1; i <= sh.size; ++i)          { sh[sh.size - i]  = norm * fc(bandwidth - i, M, 0);            }
            s = d * sh;
            for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(0, M, i) = s[i];                                    }
            
            // case f_{0,M}
            for (i = 1; i <= sh.size; ++i)          { sh[sh.size - i] = norm * fc(bandwidth - i, 0, M);             }
            if  (M & 1) { s = d * (sh * -1);} else  { s = d * sh;                                                   }
            for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(M, 0, i) = s[i];                                    }
            
            // case f_{-M,0}
            flipud(d);
            
            for (i = 1; i <= sh.size; ++i)          { sh[sh.size - i] = norm * fc(bandwidth - i, -M, 0);            }
            if (M & 1)
            {
                for (i = 0; i < sh.size; i += 2)    { sh[i] *= -1;                                                  }
            }
            else
            {
                for (i = 1; i < sh.size; i += 2)    { sh[i] *= -1;                                                  }
            }
            s = d * sh;
            for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(0, 2 * bandwidth - M, i) = s[i];                    }
            
            // case f_{0,-M}
            for (i = 1; i <= sh.size; ++i)          { sh[sh.size - i] = norm * fc(bandwidth - i, 0, -M);            }
            for (i = 1; i < sh.size; i += 2)        { sh[i] *= -1;                                                  }
            s = d * sh;
            for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(2 * bandwidth - M, 0, i) = s[i];                    }
            
            // get new wigner matrix
            d = DWT::wigner_d_matrix(bandwidth, M, M) * -1;
            d.transpose();
            
            // case f_{M,M}
            for (i = 1; i <= sh.size; ++i)          { sh[sh.size - i] = norm * fc(bandwidth - i, M, M);             }
            s = d * sh;
            for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(M, M, i) = s[i];                                    }
            
            // case f_{-M,-M}
            for (i = 1; i <= sh.size; ++i)          { sh[sh.size - i] = norm * fc(bandwidth - i, -M, -M);           }
            s = d * sh;
            for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(2 * bandwidth - M, 2 * bandwidth - M, i) = s[i];    }
            
            // Modify dw for the last two cases. flip matrix from left to right and negate every
            // second row with odd row indices.
            flipud_ne2ndocol(d);
            
            // An little arithmetic error is occuring in the following calculation... I do not exactly know why...
            // case f_{M,-M}
            for (i = 1; i <= sh.size; ++i)          { sh[sh.size - i] = norm * fc(bandwidth - i, M, -M);            }
            s = d * sh;
            for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(2 * bandwidth - M, M, i) = s[i];                    }
            
            // case f_{-M,M}
            for (i = 1; i <= sh.size; ++i)          { sh[sh.size - i] = norm * fc(bandwidth - i, -M, M);            }
            s = d * sh;
            for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(M, 2 * bandwidth - M, i) = s[i];                    }
        }
        
        // Fused two loops per hand
        //
        // for (M = 1; M < bandwidth; ++M)
        //     for (Mp = 1; Mp < M; ++Mp)
        //
        // which now is equivalent to the following loop
        #pragma omp for private(i, MMp, M, Mp, d, s, sh) schedule(dynamic) nowait
        for (MMp = 0; MMp < (bandwidth - 2) * (bandwidth - 1) / 2; ++MMp)
        {
            // reconstructing indices of the two nested for loops
            int i = MMp / (bandwidth - 1) + 1;
            int j = MMp % (bandwidth - 1) + 1;
            
            // get M and M'
            M  = j > i ? bandwidth - i : i + 1;
            Mp = j > i ? bandwidth - j : j    ;
            
            // get new wigner d-matrix
            d  = DWT::wigner_d_matrix(bandwidth, M, Mp);
            d.transpose();
            sh = vector< complex< double > >(d.cols, vec_type::COLUMN);
            
            // case f_{M,Mp}
            for (i = 1; i <= sh.size; ++i)          { sh[sh.size-i] = norm * fc(bandwidth - i, M, Mp);                  }
            sh *= -1;
            s  = d * sh;
            for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(Mp, M, i) = s[i];                                       }
            
            // case f_{Mp,M}
            for (i = 1; i <= sh.size; ++i)          { sh[sh.size - i] = norm * fc(bandwidth - i, Mp, M);                }
            if  (!((M - Mp) & 1))                   { sh *= -1;                                                         }
            s = d * sh;
            for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(M, Mp, i) = s[i];                                       }
            
            // case f_{-M,-Mp}
            for (i = 1; i <= sh.size; ++i)          { sh[sh.size - i] = norm * fc(bandwidth - i, -M, -Mp);              }
            if  (!((M - Mp) & 1))                   { sh *= -1;                                                         }
            s = d * sh;
            for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(2 * bandwidth - Mp, 2 * bandwidth - M, i) = s[i];       }
            
            // case f_{-Mp,-M}
            for (i = 1; i <= sh.size; ++i)          { sh[sh.size - i] = norm * fc(bandwidth - i, -Mp, -M);              }
            sh *= -1;
            s  = d * sh;
            for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(2 * bandwidth - M, 2 * bandwidth - Mp, i) = s[i];       }
            
            // modify wigner d-matrix for next four cases. This just works because the weight
            // function is also symmetric like the wigner-d matrix. flip up-dow the d
            // matrix and negate every second column with even row index.
            flipud_ne2ndecol(d);
            
            // case f_{Mp,-M}
            for (i = 1; i <= sh.size; ++i)          { sh[sh.size - i] = norm * fc(bandwidth - i, Mp, -M);               }
            s = d * sh;
            for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(2 * bandwidth - M, Mp, i) = s[i];                       }
            
            // case f_{M,-Mp}
            for (i = 1; i <= sh.size; ++i)          { sh[sh.size - i] = norm * fc(bandwidth - i, M, -Mp);               }
            s = d * sh;
            for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(2 * bandwidth - Mp, M, i) = s[i];                       }
            
            // alter signs
            if ((M - Mp) & 1)
            {
                for (i = 0; i < d.rows * d.cols; ++i)
                {
                    access::rw(d.mem[i]) *= -1;
                }
            }
            
            // case f_{-Mp,M}
            for (i = 1; i <= sh.size; ++i)          { sh[sh.size - i] = norm * fc(bandwidth - i, -Mp, M);               }
            s = d * sh;
            for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(M, 2 * bandwidth - Mp, i) = s[i];                       }
            
            // case f_{-M,Mp}
            for (i = 1; i <= sh.size; ++i)          { sh[sh.size - i] = norm * fc(bandwidth - i, -M, Mp);               }
            s = d * sh;
            for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(Mp, 2 * bandwidth - M, i) = s[i];                       }
        }
    }
    
    /*****************************************************************
     ** IFFT2 transform layers of input sample grid for fixed k     **
     *****************************************************************/
    synthesis.layer_wise_IDFT2(complex< double > (1. / (4. * bandwidth * bandwidth), 0));
}

UZLMATH_NAMESPACE_END

#endif
