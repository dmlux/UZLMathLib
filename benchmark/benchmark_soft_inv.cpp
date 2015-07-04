//
//  benchmark_soft_inv.cpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 14.06.15.
//
//

#ifndef benchmark_soft_for_cpp
#define benchmark_soft_for_cpp

#include <uzlmath>
#include <stdio.h>
#include <thread>

#define MAX_BW 140  // Maximal bandwidth
#define LOOP_R 10   // ISOFT runs per bandwidth

using namespace uzlmath;
using namespace FourierTransforms;

// Main method
int main(int argc, const char** argv)
{
    // space for the time values of all LOOP_R runs
    double times[MAX_BW - 1];
    
    // write to file
    char fileName[] = "benchmark_soft_inv.txt";
    FILE* fp  = fopen(fileName, "w");
    FILE* fp2 = fopen("soft_inverse.dat", "w");
    
    // print some information
    printf(     "+-----------------------------------------------------------------------------------------+\n");
    printf(     "|                                 ISOFT FORWARD BENCHMARK                                 |\n");
    printf(     "+-----------------------------------------------------------------------------------------+\n");
    printf(     "| FROM BANDWIDTH 2 TO %i WITH %i LOOP RUNS PER BANDWIDTH\n", MAX_BW, LOOP_R);
    
    // write output to file "benchmark_soft_inv.txt"
    fprintf(fp, "+-----------------------------------------------------------------------------------------+\n");
    fprintf(fp, "|                                 ISOFT FORWARD BENCHMARK                                 |\n");
    fprintf(fp, "+-----------------------------------------------------------------------------------------+\n");
    fprintf(fp, "| FROM BANDWIDTH 2 TO %i WITH %i LOOP RUNS PER BANDWIDTH\n", MAX_BW, LOOP_R);
    
    
    printf(     "+=====+===========+===================================+===================================+\n");
    printf(     "|  bw | average   | fastest run (dif. to avg / %%dif)  | slowest run (dif. to avg / %%dif)  |\n");
    printf(     "+=====+===========+===================================+===================================+\n");
    
    // write output to file "benchmark_soft_inv.txt"
    fprintf(fp, "+=====+===========+===================================+===================================+\n");
    fprintf(fp, "|  bw | average   | fastest run (dif. to avg / %%dif)  | slowest run (dif. to avg / %%dif)  |\n");
    fprintf(fp, "+=====+===========+===================================+===================================+\n");
    
    fprintf(fp2, "bandwidth\truntime\n");
    
    // loop over all bandwidth up to MAX_BW
    for (unsigned int bandwidth = 2; bandwidth <= MAX_BW; ++bandwidth)
    {
        // reset times value
        times[bandwidth - 2] = 0;
        
        // create a grid to fill with values
        grid3D< complex< double > > sample(2 * bandwidth);
        
        // create indices
        unsigned int i;
        
        // creating fourier coefficients container
        SOFTFourierCoefficients coef(bandwidth);
        
        // generate random coefficients between -1 and 1
        rand(coef, -1, 1);
        
        // min and max exec tiems
        double min, max;
        
        for (i = 0; i < LOOP_R; ++i)
        {
            // perform inverse SOFT transform
            // and stop time
            stopwatch sw = stopwatch::tic();
            ISOFT(coef, sample);
            double time  = sw.toc();
            
            // add to sum of time for current bandwidth
            times[bandwidth - 2] += time;
            
            // if first run then set min and max to the
            // first runtime
            if (i == 0)
            {
                min = time;
                max = time;
            }
            // store the current value if it is smaller than
            // the stored minimum (then it must be the new
            // minimum), or if it is bigger than max (then it
            // must be the new maximum).
            else
            {
                min = (time < min ? time : min);
                max = (time > max ? time : max);
            }
        }
        
        // get average execution time
        double avg = times[bandwidth - 2] / LOOP_R;
        
        // get fastest and slowest run-ratios
        double min_ratio = (avg - min)/avg * 100;
        double max_ratio = (max - avg)/avg * 100;
        
        // format ratios to string for right alignment
        char min_rat[7], max_rat[7];
        snprintf(min_rat, 7, "%3.3f", min_ratio);
        snprintf(max_rat, 7, "%3.3f", max_ratio);
        
        // print information
        printf("| %3i |", bandwidth);                                                // bandwidth
        printf(" %2.6fs |", avg);                                                    // average runtime for given bandwidth
        printf(" %2.6fs (-%2.6fs / -%6s%%) |",   min, (avg - min), min_rat);         // fastest run and its difference to average
        printf(" %2.6fs (+%2.6fs / +%6s%%) |\n", max, (max - avg), max_rat);         // slowest run and its difference to average
        
        // write to file
        fprintf(fp, "| %3i |", bandwidth);                                           // bandwidth
        fprintf(fp, " %2.6fs |", avg);                                               // average runtime for given bandwidth
        fprintf(fp, " %2.6fs (-%2.6fs / -%6s%%) |",   min, (avg - min), min_rat);    // fastest run and its difference to average
        fprintf(fp, " %2.6fs (+%2.6fs / +%6s%%) |\n", max, (max - avg), max_rat);    // slowest run and its difference to average
        
        fprintf(fp2, "%d\t\t%f\n", bandwidth, avg);
    }
    printf(     "+=====+===========+===================================+===================================+\n");
    fprintf(fp, "+=====+===========+===================================+===================================+\n");
    
    // close files
    fclose(fp);
    fclose(fp2);
}

#endif
