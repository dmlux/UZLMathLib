//
//  benchmark_dwt_accuracy.cpp
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
    if (argc < 5)
    {
        printf("usage: ./benchmark_dwt_for_accuracy <B> <M> <M'> <RUNS>\n");
        return 1;
    }
    
    int B    = atoi(*(argv + 1));
    int M    = atoi(*(argv + 2));
    int Mp   = atoi(*(argv + 3));
    int runs = atoi(*(argv + 4));
    
    printf("+-----------------------------------------------------------------------+\n");
    printf("|                         BENCHMARK DWT ACCURACY                        |\n");
    printf("+-----------------------------------------------------------------------+\n");
    
    printf("+------+--------------+--------------------+\n");
    printf("|  BW  | M=%d, M'=%d    | %d iterations |\n", M, Mp, runs);
    printf("+------+--------------+--------------------+\n");
    
    // run all power of 2 bandwidths
    for (int bw = 2; bw <= B; bw *= 2)
    {
        
        // create weights for the given bandwidth
        vector< double > weights = DWT::quadrature_weights(bw);
    
        // create weigthed wigner matrix and wigner matrix
        matrix< double > dw = DWT::weighted_wigner_d_matrix(bw, M, Mp, weights);
        matrix< double > dt = DWT::wigner_d_matrix(bw, M, Mp);
        dt.transpose();
        
        // relative and absolute errors
        double relative = 0;
        double absolute = 0;
        
        // run several runs and take average values
        for (int i = 0; i < runs; ++i)
        {
            
            // create random coefficients
            vector< complex< double > > fh(dt.cols, vec_type::COLUMN);
            rand(fh, -1, 1);
            
            // inverse DWT
            vector< complex< double > > s = dt * fh;
            
            // forward DWT
            vector< complex< double > > gh = dw * s;
            
            // get difference vector
            vector< complex< double > > dif = gh - fh;
            
            // get maximum
            double max_val = dif[0].abs();
            int max_idx = 0;
            
            for (int j = 0; j < dif.size; ++j)
            {
                if (dif[j].abs() > max_val)
                {
                    max_val = dif[j].abs();
                    max_idx = j;
                }
            }
            
            relative += max_val / fh[max_idx].abs();
            absolute += max_val;
        }
        
        relative /= runs;
        absolute /= runs;
        
        printf("| %4d | %e | absolute error     |\n", bw, absolute);
        printf("|      | %e | relative error     |\n", relative);
        printf("+------+--------------+--------------------+\n");
    }
    
    return 0;
}