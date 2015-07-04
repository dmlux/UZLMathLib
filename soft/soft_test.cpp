//
//  SOFT.cpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 25.05.15.
//
//

#ifndef SOFT_cpp
#define SOFT_cpp

#include <uzlmath>
#include <stdio.h>
#include <chrono>
#include <thread>
#include <fftw3.h>

#if _OPENMP
    #include <omp.h>
#endif

using namespace uzlmath;

complex< double > f(const double& alpha, const double& beta, const double& gamma)
{
    return -7.31 * wigner::wigner_D_l2normalized(2,  0,  -1, -alpha, -beta, -gamma)
    + complex<double>(1., -9.731) * wigner::wigner_D_l2normalized(1,  0,  0, -alpha, -beta, -gamma);
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
    FourierTransforms::SOFT(sample, coef);
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
    FourierTransforms::ISOFT(coef, grid_rec);
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
    FourierTransforms::ISOFT(coef, sample);
    double time2 = sw.toc();
    
    // perform forward SOFT transform
    sw = stopwatch::tic();
    FourierTransforms::SOFT(sample, rec_coef);
    double time = sw.toc();
    
    // print Fourier coefficients
    // save outstream flags
    if (show_coefs)
    {
        printf("** Fourier coefficients:\n");
    }
    
    bool equal = true;
    double epsilon = 1e-11;
    
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
                
                if ( show_coefs )
                {
                    // Here the coefficients are printed out on the console
                    printf("l=%4d, M=%4d, M'=%4d: %.4f%s%.4f\n", m,  n, k, coef(m,n,k).re, (coef(m,n,k).im >= 0 ? "+" : ""), coef(m,n,k).im);
                    printf("%.16f\n", fabs(coef(m,n,k).re - rec_coef(m,n,k).re));
                    
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
    
    // store coefficients on disk
//    obj2file(rec_coef, "fc.txt");
    
    printf("Bandbreite:     %d\n", bandwidth);
    printf("SOFT:           %.6fs\n", time);
    printf("ISOFT:          %.6fs\n", time2);
    printf("Correct result: %s\n", (equal ? "Yes" : "No"));
}

int main(int argc, const char ** argv)
{
    //createGridSOFT(10);
//    for_back_file("/Users/dlux/Desktop/soft_files/grid_128_samp.dat", 128, false);
//    for_back_file("/Users/dlux/Desktop/soft_files/test_series/grid_3_test.dat", 3, true);
//    for_back(128, false);
    
    int n = 1024 * 2;
    matrix< complex< double > > A(n), B(n), C(n);
    
    rand(A, -1, 1);
    B = A;
    
    stopwatch sw = stopwatch::tic();
    FourierTransforms::DFT2(B);
    double time1 = sw.toc();
    
    // create a FFTW plan and executing an FFT2 by use of 1D ffts
    // First part is the none strided version... (columns)
    memcpy(C.memptr(), A.memptr(), A.n_cols() * A.n_rows() * sizeof(complex< double >));
    
    sw = stopwatch::tic();
    fftw_complex* mat = reinterpret_cast< fftw_complex* >(C.memptr());
    fftw_plan* cols = new fftw_plan[C.n_cols()], * rows = new fftw_plan[C.n_rows()];
    
    // first transform the cols of the matrix
    for (int i = 0; i < C.n_cols(); ++i)
    {
        cols[i] = fftw_plan_dft_1d(C.n_rows(), mat + i * C.n_rows(), mat + i * C.n_rows(), FFTW_FORWARD, FFTW_ESTIMATE);
    }
    
    // second transform the rows of the matrix
    for (int i = 0; i < C.n_rows(); ++i)
    {
        int N = C.n_cols();
        int istride = C.n_rows();
        rows[i] = fftw_plan_many_dft(1, &N, 1, mat + i, NULL, istride, 1, mat + i, NULL, istride, 1, FFTW_FORWARD, FFTW_ESTIMATE);
    }
    
    // execute column FFTs
    #pragma omp parallel for shared(cols, C) schedule(dynamic)
    for (int i = 0; i < C.n_cols(); ++i)
    {
        fftw_execute(cols[i]);
    }
    
    // execute row FFTs
    #pragma omp parallel for shared(rows, C) schedule(dynamic)
    for (int i = 0; i < C.n_rows(); ++i)
    {
        fftw_execute(rows[i]);
    }
    
    // destroy plans and cleanup
    for (int i = 0; i < C.n_rows(); ++i)
    {
        fftw_destroy_plan(rows[i]);
    }
    for (int i = 0; i < C.n_cols(); ++i)
    {
        fftw_destroy_plan(cols[i]);
    }
    delete [] rows;
    delete [] cols;
    fftw_cleanup();
    double time2 = sw.toc();
    
    bool equal = true;
    for (int i = 0; i < A.n_rows() * A.n_cols(); ++i)
    {
        if (fabs(B.memptr()[i].re - C.memptr()[i].re) > 1e-10 || fabs(B.memptr()[i].im - C.memptr()[i].im) > 1e-10)
        {
            equal = false;
            break;
        }
    }
    
    std::cout << "time1: " << std::fixed << std::setprecision(6) << time1 << "s" << std::endl << "time2: " << time2 << "s" << std::endl;
    std::cout << "Equal? " << (equal ? "Yes" : "No") << std::endl;
    
//    std::cout << "B = " << B << std::endl;
//    std::cout << "C = " << C << std::endl;
    
    return 0;
}

#endif
