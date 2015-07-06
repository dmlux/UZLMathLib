//
//  fn_fourier_transforms.cpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 18.06.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_fn_fourier_transforms_cpp
#define uzlmath_fn_fourier_transforms_cpp

#include <uzlmath>

// include OpenMP if compiler supports it
#if _OPENMP
    #include <omp.h>
#endif

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
auto SOFT(grid3D< complex< double > > sample, SOFTFourierCoefficients& fc) -> void
{
    /*****************************************************************
     ** Check parameters                                            **
     *****************************************************************/
    // Check if the grid has same size in each dimension
    if (sample.n_rows() != sample.n_cols() || sample.n_rows() != sample.n_lays())
    {
        printf("** uzlmath error: all SOFT sample grid dimensions should be equal. **\n");
        return;
    }
    
    // Check if grid has odd dimensions
    if (sample.n_rows() & 1)
    {
        printf("** uzlmath error: SOFT sample grid dimensions are not even. **\n");
        return;
    }
    
    // Extract bandwidth
    unsigned int bandwidth = static_cast< unsigned int >(sample.n_cols() / 2);
    
    // Check if Fourier coefficients container dimension matches sample dimension
    if (bandwidth != fc.bandwidth)
    {
        printf("** uzlmath error: SOFT Fourier coefficients container bandwidth does not match to sample grid bandwidth. **\n");
        return;
    }
    
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
    double norm              = M_PI / (2 * bandwidth * bandwidth);
    
    // defining needed indices
    unsigned int MMp, i, M, Mp;
    
    // DWT for M = 0, M' = 0
    for (i = 0; i < 2 * bandwidth; ++i)     { s[i] = sample(0, 0, i);                               }
    vector< complex< double > > sh = dw * s;
    for (i = 1; i <= sh.n_elements(); ++i)  { fc(bandwidth-i, 0, 0) = norm * sh[sh.n_elements()-i]; }
    
    /*****************************************************************
     ** Iterate over all combinations of M and M'                   **
     *****************************************************************/
    #pragma omp parallel for default(shared) private(i, M, dw, sh) firstprivate(s) schedule(dynamic)
    for (M = 1; M < bandwidth; ++M)
    {
        dw = DWT::weighted_wigner_d_matrix(bandwidth, M, 0, weights) * -1;
        
        /*****************************************************************
         ** Make use of symmetries                                      **
         *****************************************************************/
        // case f_{M,0}
        for (i = 0; i < 2 * bandwidth; ++i)     { s[i] = sample(0, M, i);                                 }
        sh = dw * s;
        for (i = 1; i <= sh.n_elements(); ++i)  { fc(bandwidth-i, M, 0) = norm * sh[sh.n_elements()-i];   }
        
        // case f_{0,M}
        for (i = 0; i < 2 * bandwidth; ++i)     { s[i] = sample(M, 0, i);                                 }
        sh = dw * s;
        if  (M & 1)                             { sh *= -1;                                               }
        for (i = 1; i <= sh.n_elements(); ++i)  { fc(bandwidth-i, 0, M) = norm * sh[sh.n_elements()-i];   }
        
        // case f_{-M,0}
        fliplr(dw);
        for (i = 0; i < 2 * bandwidth; ++i)     { s[i] = sample(0, 2 * bandwidth - M, i);                 }
        sh = dw * s;
        if (M & 1)
        {
            for (i = 0; i < sh.n_elements(); i += 2) { sh[i] *= -1;                                       }
        }
        else
        {
            for (i = 1; i < sh.n_elements(); i += 2) { sh[i] *= -1;                                       }
        }
        for (i = 1; i <= sh.n_elements(); ++i)  { fc(bandwidth-i, -M, 0) = norm * sh[sh.n_elements()-i];  }
        
        // case f_{0,-M}
        for (i = 0; i < 2 * bandwidth; ++i)     { s[i] = sample(2 * bandwidth - M, 0, i);                 }
        sh = dw * s;
        for (i = 1; i < sh.n_elements(); i +=2 ){ sh[i] *= -1;                                            }
        for (i = 1; i <= sh.n_elements(); ++i)  { fc(bandwidth-i, 0, -M) = norm * sh[sh.n_elements()-i];  }
        
        // get new wigner matrix
        dw = DWT::weighted_wigner_d_matrix(bandwidth, M, M, weights) * -1;
        
        // case f_{M, M}
        for (i = 0; i < 2 * bandwidth; ++i)     { s[i] = sample(M, M, i);                                 }
        sh = dw * s;
        for (i = 1; i <= sh.n_elements(); ++i)  { fc(bandwidth-i, M, M) = norm * sh[sh.n_elements()-i];   }
        
        // case f_{-M, -M}
        for (i = 0; i < 2 * bandwidth; ++i)     { s[i] = sample(2 * bandwidth - M, 2 * bandwidth - M, i); }
        sh = dw * s;
        for (i = 1; i <= sh.n_elements(); ++i)  { fc(bandwidth-i, -M, -M) = norm * sh[sh.n_elements()-i]; }
        
        // Modify dw for the last two cases. flip matrix from left to right and negates sign of
        // every second row with odd row index.
        fliplr_ne2ndorow(dw);
        
        // An little arithmetic error is occuring in the following calculation
        // case f_{M, -M}
        for (i = 0; i < 2 * bandwidth; ++i)     { s[i] = sample(2 * bandwidth - M, M, i);                 }
        sh = dw * s;
        for (i = 1; i <= sh.n_elements(); ++i)  { fc(bandwidth-i, M, -M) = norm * sh[sh.n_elements()-i];  }
        
        // case f_{-M, M}
        for (i = 0; i < 2 * bandwidth; ++i)     { s[i] = sample(M, 2 * bandwidth - M, i);                 }
        sh = dw * s;
        for (i = 1; i <= sh.n_elements(); ++i)  { fc(bandwidth-i, -M, M) = norm * sh[sh.n_elements()-i];  }
    }
    
    // Fused two loops per hand
    //
    // for (M = 1; M < bandwidth; ++M)
    //     for (Mp = 1; Mp < M; ++Mp)
    //
    // which now is equivalent to the following loop
    #pragma omp parallel for default(shared) private(i, MMp, M, Mp, dw, sh) firstprivate(s) schedule(dynamic)
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
        for (i = 0; i < 2 * bandwidth; ++i)     { s[i] = sample(Mp, M, i);                                          }
        sh  = dw * s;
        sh *= -1;
        for (i = 1; i <= sh.n_elements(); ++i)  { fc(bandwidth-i, M, Mp) = norm * sh[sh.n_elements()-i];            }
        
        // case f_{Mp, M}
        for (i = 0; i < 2 * bandwidth; ++i)     { s[i] = sample(M, Mp, i);                                          }
        sh = dw * s;
        if  (!((M - Mp) & 1))                   { sh *= -1;                                                         }
        for (i = 1; i <= sh.n_elements(); ++i)  { fc(bandwidth-i, Mp, M) = norm * sh[sh.n_elements()-i];            }
        
        // case f_{-M, -Mp}
        for (i = 0; i < 2 * bandwidth; ++i)     { s[i] = sample(2 * bandwidth - Mp, 2 * bandwidth - M, i);          }
        sh = dw * s;
        if  (!((M - Mp) & 1))                   { sh *= -1;                                                         }
        for (i = 1; i <= sh.n_elements(); ++i)  { fc(bandwidth-i, -M, -Mp) = norm * sh[sh.n_elements()-i];          }
        
        // case f_{-Mp, -M}
        for (i = 0; i < 2 * bandwidth; ++i)     { s[i] = sample(2 * bandwidth - M, 2 * bandwidth - Mp, i);          }
        sh  = dw * s;
        sh *= -1;
        for (i = 1; i <= sh.n_elements(); ++i)  { fc(bandwidth-i, -Mp, -M) = norm * sh[sh.n_elements()-i];          }
        
        // modify wigner d-matrix for next four cases. This just works because the weight
        // function is also symmetric like the wigner-d matrix. flip left-right the dw
        // matrix and negate each even value with even row index.
        fliplr_ne2nderow(dw);
        
        // case f_{Mp, -M}
        for (i = 0; i < 2 * bandwidth; ++i)     { s[i] = sample(2 * bandwidth - M, Mp, i);                          }
        sh = dw * s;
        for (i = 1; i <= sh.n_elements(); ++i)  { fc(bandwidth-i, Mp, -M) = norm * sh[sh.n_elements()-i];           }
        
        // case f_{M, -Mp}
        for (i = 0; i < 2 * bandwidth; ++i)     { s[i] = sample(2 * bandwidth - Mp, M, i);                          }
        sh = dw * s;
        for (i = 1; i <= sh.n_elements(); ++i)  { fc(bandwidth-i, M, -Mp) = norm * sh[sh.n_elements()-i];           }
        
        // alter signs
        if ((M - Mp) & 1)
        {
            double* dw_mem = dw.memptr();
            for (i = 0; i < dw.n_rows() * dw.n_cols(); ++i)
            {
                dw_mem[i] *= -1;
            }
        }
        
        // case f_{-Mp, M}
        for (i = 0; i < 2 * bandwidth; ++i)     { s[i] = sample(M, 2 * bandwidth - Mp, i);                          }
        sh = dw * s;
        for (i = 1; i <= sh.n_elements(); ++i)  { fc(bandwidth-i, -Mp, M) = norm * sh[sh.n_elements()-i];           }
        
        // case f_{-M, Mp}
        for (i = 0; i < 2 * bandwidth; ++i)     { s[i] = sample(Mp, 2 * bandwidth - M, i);                          }
        sh = dw * s;
        for (i = 1; i <= sh.n_elements(); ++i)  { fc(bandwidth-i, -M, Mp) = norm * sh[sh.n_elements()-i];           }
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
auto ISOFT(const SOFTFourierCoefficients& fc, grid3D< complex< double > >& synthesis) -> void
{
    /*****************************************************************
     ** Check parameters                                            **
     *****************************************************************/
    // Check if the grid has same size in each dimension
    if (synthesis.n_rows() != synthesis.n_cols() || synthesis.n_rows() != synthesis.n_lays())
    {
        printf("** uzlmath error: all ISOFT synthesis grid dimensions should be equal. **\n");
        return;
    }
    
    // Check if grid has odd dimensions
    if (synthesis.n_rows() & 1)
    {
        printf("** uzlmath error: ISOFT synthesis grid dimensions are not even. **\n");
        return;
    }
    
    // Extract bandwidth
    unsigned int bandwidth = static_cast< unsigned int >(synthesis.n_cols() / 2);
    
    // Check if Fourier coefficients container dimension matches sample dimension
    if (bandwidth != fc.bandwidth)
    {
        printf("** uzlmath error: ISOFT Fourier coefficients container bandwidth does not match to synthesis grid bandwidth. **\n");
        return;
    }
    
    /*****************************************************************
     ** M = 0, M' = 0                                               **
     *****************************************************************/
    matrix< double > d = DWT::wigner_d_matrix(bandwidth, 0, 0) * -1;
    vector< complex< double > > sh(d.n_rows(), vec_type::COLUMN);
    
    d.transpose();
    
    // defining norm factor
    double norm = (2 * bandwidth * bandwidth) / M_PI;
    
    // defining needed indices
    unsigned MMp, i, M, Mp;
    
    // inverse DWT for M = 0, M' = 0
    for (i = 1; i <= sh.n_elements(); ++i)  { sh[sh.n_elements()-i] = norm * fc(bandwidth - i, 0, 0); }
    vector< complex< double > > s = d * sh;
    for (i = 0; i < 2 * bandwidth; ++i) { synthesis(0, 0, i) = s[i];                                  }
    
    /*****************************************************************
     ** Iterate over all combinations of M and M'                   **
     *****************************************************************/
    #pragma omp parallel for private(i, M, d, s, sh) schedule(dynamic)
    for (M = 1; M < bandwidth; ++M)
    {
        d  = DWT::wigner_d_matrix(bandwidth, M, 0) * -1;
        sh = vector< complex< double > >(d.n_rows(), vec_type::COLUMN);
        d.transpose();
        
        /*****************************************************************
         ** Make use of symmetries                                      **
         *****************************************************************/
        // case f_{M,0}
        for (i = 1; i <= sh.n_elements(); ++i)  { sh[sh.n_elements()-i]  = norm * fc(bandwidth-i, M, 0);  }
        s = d * sh;
        for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(0, M, i) = s[i];                              }
        
        // case f_{0,M}
        for (i = 1; i <= sh.n_elements(); ++i)  { sh[sh.n_elements()-i] = norm * fc(bandwidth-i, 0, M);   }
        if  (M & 1) { s = d * (sh * -1);} else  { s = d * sh;                                             }
        for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(M, 0, i) = s[i];                              }
        
        // case f_{-M,0}
        flipud(d);
        
        for (i = 1; i <= sh.n_elements(); ++i)  { sh[sh.n_elements()-i] = norm * fc(bandwidth-i, -M, 0);  }
        if (M & 1)
        {
            for (i = 0; i < sh.n_elements(); i += 2) { sh[i] *= -1;                                       }
        }
        else
        {
            for (i = 1; i < sh.n_elements(); i += 2) { sh[i] *= -1;                                       }
        }
        s = d * sh;
        for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(0, 2 * bandwidth - M, i) = s[i];              }
        
        // case f_{0,-M}
        for (i = 1; i <= sh.n_elements(); ++i)  { sh[sh.n_elements()-i] = norm * fc(bandwidth-i, 0, -M);  }
        for (i = 1; i < sh.n_elements(); i += 2){ sh[i] *= -1;                                            }
        s = d * sh;
        for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(2 * bandwidth - M, 0, i) = s[i];              }
        
        // get new wigner matrix
        d = DWT::wigner_d_matrix(bandwidth, M, M) * -1;
        d.transpose();
        
        // case f_{M,M}
        for (i = 1; i <= sh.n_elements(); ++i)  { sh[sh.n_elements()-i] = norm * fc(bandwidth-i, M, M);   }
        s = d * sh;
        for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(M, M, i) = s[i];                              }
        
        // case f_{-M,-M}
        for (i = 1; i <= sh.n_elements(); ++i)  { sh[sh.n_elements()-i] = norm * fc(bandwidth-i, -M, -M); }
        s = d * sh;
        for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(2 * bandwidth - M, 2 * bandwidth - M, i) = s[i]; }
        
        // Modify dw for the last two cases. flip matrix from left to right and negate every
        // second row with odd row index.
        flipud_ne2ndocol(d);
        
        // An little arithmetic error is occuring in the following calculation
        // case f_{M,-M}
        for (i = 1; i <= sh.n_elements(); ++i)  { sh[sh.n_elements()-i] = norm * fc(bandwidth-i, M, -M);  }
        s = d * sh;
        for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(2 * bandwidth - M, M, i) = s[i];              }
        
        // case f_{-M,M}
        for (i = 1; i <= sh.n_elements(); ++i)  { sh[sh.n_elements()-i] = norm * fc(bandwidth-i, -M, M);  }
        s = d * sh;
        for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(M, 2 * bandwidth - M, i) = s[i];              }
    }
    
    // Fused two loops per hand
    //
    // for (M = 1; M < bandwidth; ++M)
    //     for (Mp = 1; Mp < M; ++Mp)
    //
    // which now is equivalent to the following loop
    #pragma omp parallel for private(i, MMp, M, Mp, d, s, sh) schedule(dynamic)
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
        sh = vector< complex< double > >(d.n_cols(), vec_type::COLUMN);
        
        // case f_{M,Mp}
        for (i = 1; i <= sh.n_elements(); ++i)  { sh[sh.n_elements()-i] = norm * fc(bandwidth-i, M, Mp);            }
        sh *= -1;
        s  = d * sh;
        for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(Mp, M, i) = s[i];                                       }
        
        // case f_{Mp,M}
        for (i = 1; i <= sh.n_elements(); ++i)  { sh[sh.n_elements()-i] = norm * fc(bandwidth-i, Mp, M);            }
        if  (!((M - Mp) & 1))                   { sh *= -1;                                                         }
        s = d * sh;
        for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(M, Mp, i) = s[i];                                       }
        
        // case f_{-M,-Mp}
        for (i = 1; i <= sh.n_elements(); ++i)  { sh[sh.n_elements()-i] = norm * fc(bandwidth-i, -M, -Mp);          }
        if  (!((M - Mp) & 1))                   { sh *= -1;                                                         }
        s = d * sh;
        for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(2 * bandwidth - Mp, 2 * bandwidth - M, i) = s[i];       }
        
        // case f_{-Mp,-M}
        for (i = 1; i <= sh.n_elements(); ++i)  { sh[sh.n_elements()-i] = norm * fc(bandwidth-i, -Mp, -M);          }
        sh *= -1;
        s  = d * sh;
        for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(2 * bandwidth - M, 2 * bandwidth - Mp, i) = s[i];       }
        
        // modify wigner d-matrix for next four cases. This just works because the weight
        // function is also symmetric like the wigner-d matrix. flip up-dow the d
        // matrix and negate every second column with even row index.
        flipud_ne2ndecol(d);
        
        // case f_{Mp,-M}
        for (i = 1; i <= sh.n_elements(); ++i)  { sh[sh.n_elements()-i] = norm * fc(bandwidth-i, Mp, -M);           }
        s = d * sh;
        for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(2 * bandwidth - M, Mp, i) = s[i];                       }
        
        // case f_{M,-Mp}
        for (i = 1; i <= sh.n_elements(); ++i)  { sh[sh.n_elements()-i] = norm * fc(bandwidth-i, M, -Mp);           }
        s = d * sh;
        for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(2 * bandwidth - Mp, M, i) = s[i];                       }
        
        // alter signs
        if ((M - Mp) & 1)
        {
            double* d_mem = d.memptr();
            for (i = 0; i < d.n_rows() * d.n_cols(); ++i)
            {
                d_mem[i] *= -1;
            }
        }
        
        // case f_{-Mp,M}
        for (i = 1; i <= sh.n_elements(); ++i)  { sh[sh.n_elements()-i] = norm * fc(bandwidth-i, -Mp, M);           }
        s = d * sh;
        for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(M, 2 * bandwidth - Mp, i) = s[i];                       }
        
        // case f_{-M,Mp}
        for (i = 1; i <= sh.n_elements(); ++i)  { sh[sh.n_elements()-i] = norm * fc(bandwidth-i, -M, Mp);           }
        s = d * sh;
        for (i = 0; i < 2 * bandwidth; ++i)     { synthesis(Mp, 2 * bandwidth - M, i) = s[i];                       }
    }
    
    /*****************************************************************
     ** IFFT2 transform layers of input sample grid for fixed k     **
     *****************************************************************/
    synthesis.layer_wise_IDFT2(complex< double > (1.0 /(4 * bandwidth * bandwidth), 0));
}

UZLMATH_NAMESPACE_END

#endif
