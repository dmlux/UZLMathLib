//
//  SOFT.cpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 25.05.15.
//
//

#ifndef SOFFT_cpp
#define SOFFT_cpp

#include <uzlmath>
#include <stdio.h>
#include <chrono>
#include <thread>

using namespace uzlmath;

uzlmath_deprecated
complex< double > f(const double& alpha, const double& beta, const double& gamma)
{
    return complex< double >(-7.31, 0     ) * wigner::wigner_D_l2normalized(2,  0,  -1, -alpha, -beta, -gamma)
         + complex< double >(1.   , -9.731) * wigner::wigner_D_l2normalized(1,  0,  0, -alpha, -beta, -gamma);
//        + 13 * wigner::wigner_D_l2normalized(2,  1,  0, -alpha, -beta, -gamma)
//        - 8.423 * wigner::wigner_D_l2normalized(2,  0,  2, -alpha, -beta, -gamma)
//        + complex<double>(1.8312, -9.8372) * wigner::wigner_D_l2normalized(2,  0,  0, -alpha, -beta, -gamma);
}

void createGridSOFT(unsigned int B)
{
    // create grid
    grid3D< complex< double > > grid(2 * B, 0);
    
    // filling grid with custom wigners given in f
    int cnt, l = 0, M, Mp;
    for (int k = 0; k < 2 * B; ++k)
    {
        for (int i = 0; i < 2 * B; ++i)
        {
            for (int j = 0; j < 2 * B; ++j)
            {
                // grid(i, j, k) = f((2*M_PI*j)/(2*B), (M_PI * (2.0 * k + 1.0)) / (4.0 * B), (2*M_PI*i)/(2*B));
                
                // use coefficients in numerical order
                cnt = 1;
                for (l = 0; l < B; ++l)
                {
                    for (M = -l; M <= l; ++M)
                    {
                        for (Mp = -l; Mp <= l; ++Mp)
                        {
                            grid(i, j, k) += complex< double >(cnt, cnt+1) * wigner::wigner_D_l2normalized(l, M, Mp, -(2*M_PI*j)/(2*B), -(M_PI * (2.0 * k + 1.0)) / (4.0 * B), -(2*M_PI*i)/(2*B));
                            cnt += 2;
                        }
                    }
                }
            }
        }
    }
    
    
    // write to file
    char fileName[100];
    std::snprintf(fileName, 100, "grid_%i_samp.dat", B);
    FILE* fp = fopen(fileName, "w");
    for (int i = 0; i < 2 * B; ++i)
    {
        for (int j = 0; j < 2 * B; ++j)
        {
            for (int k = 0; k < 2 * B; ++k)
            {
                fprintf(fp, "%.16f\n%.16f\n", grid(j, k, i).re, grid(j, k, i).im);
            }
        }
    }
    fclose(fp);
}

void for_back_file(const char* fileName, unsigned int bandwidth, bool show_coefs)
{
    // create a grid to fill with values
    grid3D< complex< double > > sample(2 * bandwidth);
    
    // create indices
    unsigned int m = 0, n = 0, k = 0;
    
    // read doubles from file
    stopwatch sw = stopwatch::tic();
    FILE* fp = fopen( fileName, "r" );
    for ( unsigned int i = 0 ; i < 8 * bandwidth * bandwidth * bandwidth; i++ )
    {
        double re, im;
        
        fscanf(fp, "%lf", &re);
        fscanf(fp, "%lf", &im);
        
        sample(m, k, n).re = re;
        sample(m, k, n).im = im;
        
        k++;
        if (k >= 2 * bandwidth)
        {
            k = 0;
            m++;
            if (m >= 2 * bandwidth)
            {
                m = 0;
                n++;
            }
        }
    }
    fclose ( fp ) ;
    double time = sw.toc();
    
    // creating fourier coefficients container
    SOFTFourierCoefficients coef(bandwidth);
    
    // perform forward SOFT transform
    sw = stopwatch::tic();
    FourierTransforms::DSOFT(sample, coef);
    time = sw.toc();
    
    // print Fourier coefficients
    // save outstream flags
    if (show_coefs)
    {
        printf("** Fourier coefficients:\n");
    }
    for (int m = 0; m < bandwidth; ++m)
    {
        for (int n = -m; n <= m; ++n)
        {
            for (int k = -m; k <= m; ++k)
            {
                if ( show_coefs )
                {
                    // Here the coefficients are printed out on the console
                    printf("l=%4d, M=%4d, M'=%4d: %.4f%s%.4f\n", m,  n, k, coef(m,n,k).re, (coef(m,n,k).im >= 0 ? "+" : ""), coef(m,n,k).im);
                }
            }
        }
    }
    if (show_coefs)
    {
        printf("\n");
    }
    
    // store coefficients on disk
    obj2file(coef, "fc.txt");
    
    grid3D< complex< double > > grid_rec(2 * bandwidth);
    
    sw = stopwatch::tic();
    FourierTransforms::IDSOFT(coef, grid_rec);
    double time2 = sw.toc();
    
    printf("Bandbreite:   %d\n", bandwidth);
    printf("SOFT:         %.6fs\n", time);
    printf("ISOFT:        %.6fs\n", time2);
    
    std::cout << "sample = " << sample << std::endl;
    std::cout << "reconstructed sample = " << grid_rec << std::endl;
}

void for_back(unsigned int bandwidth, bool show_coefs)
{
    // create a grid to fill with values
    grid3D< complex< double > > sample(2 * bandwidth);
    
    // creating fourier coefficients container
    SOFTFourierCoefficients coef(bandwidth);
    SOFTFourierCoefficients rec_coef(bandwidth);
    
    rand(coef, -1, 1);
    
    stopwatch sw = stopwatch::tic();
    FourierTransforms::IDSOFT(coef, sample);
    double time2 = sw.toc();
    
    // perform forward SOFT transform
    sw = stopwatch::tic();
    FourierTransforms::DSOFT(sample, rec_coef);
    double time = sw.toc();
    
//    std::cout << "coef = " << coef << std::endl;
//    std::cout << "rec_coef = " << rec_coef << std::endl;
//    std::cout << sample << std::endl;
    
    // print Fourier coefficients
    // save outstream flags
    if (show_coefs)
    {
        printf("** Fourier coefficients:\n");
    }
    
    bool equal = true;
    double epsilon = 1e-11;
    int cnt_fc = 0;
    
    for (int m = 0; m < bandwidth; ++m)
    {
        for (int n = -m; n <= m; ++n)
        {
            for (int k = -m; k <= m; ++k)
            {
//                printf("%.16f\n", fabs(coef(m,n,k).re - rec_coef(m,n,k).re));
                if (fabs(coef(m,n,k).re - rec_coef(m,n,k).re) > epsilon || fabs(coef(m,n,k).im - rec_coef(m,n,k).im) > epsilon)
                {
                    equal = false;
                }
                
                cnt_fc++;
                
                if ( show_coefs )
                {
                    // Here the coefficients are printed out on the console
                    printf("l=%4d, M=%4d, M'=%4d: %.4f%s%.4f\n", m,  n, k, coef(m,n,k).re, (coef(m,n,k).im >= 0 ? "+" : ""), coef(m,n,k).im);
                    //printf("%.16f\n", fabs(coef(m,n,k).re - rec_coef(m,n,k).re));
                    
                    if (fabs(coef(m,n,k).re - rec_coef(m,n,k).re) > epsilon || fabs(coef(m,n,k).im - rec_coef(m,n,k).im) > epsilon)
                    {
                        equal = false;
                    }
                }
            }
        }
    }
    if (show_coefs)
    {
        printf("\n");
    }
    
    printf("#coefficients:  %d\n", cnt_fc);
    
    // store coefficients on disk
//    obj2file(rec_coef, "fc.txt");
    
    printf("Bandbreite:     %d\n", bandwidth);
    printf("SOFT:           %.6fs\n", time);
    printf("ISOFT:          %.6fs\n", time2);
    printf("Correct result: %s\n", (equal ? "Yes" : "No"));
}

int main(int argc, const char ** argv)
{
    if (argc < 2)
    {
        printf("usage: ./soft_test <Bandwidth>\n");
        return 1;
    }
    
    int B = atoi(*(argv + 1));
        
    //createGridSOFT(10);
//    for_back_file("/Users/dlux/Desktop/soft_files/grid_128_samp.dat", 128, false);
//    for_back_file("/Users/dlux/Desktop/soft_files/test_series/grid_3_test.dat", 3, true);
    for_back(B, false);
    
    return 0;
}

#endif
