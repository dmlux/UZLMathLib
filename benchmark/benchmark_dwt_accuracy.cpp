//
//  benchmark_soft_for.cpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 14.06.15.
//
//

#include <uzlmath>
#include <stdio.h>

using namespace uzlmath;

int main(int argc, const char** argv)
{
    if (argc < 3)
    {
        printf("usage: ./benchmark_dwt_for_accuracy <B> <RUNS>\n");
        return 1;
    }
    
    int B    = atoi(*(argv + 1));
    int runs = atoi(*(argv + 2));
    
    // create weights for the given bandwidth
    vector< double > weights = DWT::quadrature_weights(B);
    
    printf("+-----------------------------------------------------------------------+\n");
    printf("|                         BENCHMARK DWT FORWARD                         |\n");
    printf("+-----------------------------------------------------------------------+\n");
    
    for (int M = -(B - 1); M < B; ++M)
    {
        for (int Mp = -(B - 1); Mp < B; ++Mp)
        {
            // create weigthed wigner matrix and wigner matrix
            matrix< double > dw = DWT::weighted_wigner_d_matrix(B, M, Mp, weights);
            matrix< double > d  = DWT::wigner_d_matrix(B, M, Mp);
            d.transpose();
            
            // create random coefficients
            vector< complex< double > > sh(d.cols, vec_type::COLUMN);
            rand(sh, -1, 1);
            
            vector< complex< double > > sh_tmp = sh;
            vector< complex< double > > s_tmp(2 * B, vec_type::COLUMN);
            
            // run forward and inverse transform
            for (int i = 0; i < runs; ++i)
            {
                s_tmp = d * sh_tmp;
                sh_tmp = dw * s_tmp;
            }
            
            // difference vector
            vector< complex< double > > dif = sh - sh_tmp;
            
            // get maximum value
            complex< double > max = dif[0];
            
            // get maximum in differences
            for (int i = 0; i < dif.size; ++i)
            {
                if (max.abs() > dif[i].abs())
                {
                    max = dif[i];
                }
            }
            
            printf("M: %d, M': %d, max error: %e + %ei\n", M, Mp, max.re, max.im);
        }
    }
    
    return 0;
}