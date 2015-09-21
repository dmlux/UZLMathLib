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
    
    printf("+------+--------------+----------------+\n");
    printf("|  BW  | M=%d, M'=%d    | %d iterations |\n", M, Mp, runs);
    printf("+------+--------------+----------------+\n");
    
    // run all power of 2 bandwidths
    for (int bw = 2; bw <= B; bw *= 2)
    {
        
        // create weights for the given bandwidth
        vector< double > weights = DWT::quadrature_weights(bw);
    
        // create weigthed wigner matrix and wigner matrix
        matrix< double > dw = DWT::weighted_wigner_d_matrix(bw, M, Mp, weights);
        matrix< double > dt = DWT::wigner_d_matrix(bw, M, Mp);
        dt.transpose();
        
        // absolute average error
        double abs_exact = 0;
        double abs_error = 0;
        
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
            double max_error = dif[0].abs();
            double max_exact =  fh[0].abs();
            
            // getting max inifinity norm
            for (int j = 0; j < dif.size; ++j)
            {
                // difference vector
                if (dif[j].abs() > max_error)
                {
                    max_error = dif[j].abs();
                }
                
                // exact vector
                if (fh[j].abs() > max_exact)
                {
                    max_exact =  fh[j].abs();
                }
            }
            
            // add to absolute error
            abs_error += max_error;
            abs_exact += max_exact;
        }
        
        // divide abs by the amount of runs to get average
        // abs error
        abs_error /= runs;
        abs_exact /= runs;
        
        //printf("%e, %e\n", abs_error, abs_exact);
        
        printf("| %4d | %e | absolute error |\n", bw, abs_error);
        printf("|      | %e | relative error |\n", abs_error / abs_exact);
        printf("+------+--------------+----------------+\n");
    }
    
    return 0;
}