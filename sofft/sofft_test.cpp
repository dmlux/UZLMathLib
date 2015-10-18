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

//uzlmath_deprecated
complex< double > f(const double& alpha, const double& beta, const double& gamma)
{
    return
          complex< double >(-7.31, 0)       * Wigner::wigner_D_l2normalized(2,  0, -1, -alpha, -beta, -gamma)
        + complex< double >(1, -9.731)      * Wigner::wigner_D_l2normalized(1,  0,  0, -alpha, -beta, -gamma)
        + complex< double >(13, 0)          * Wigner::wigner_D_l2normalized(2,  1,  0, -alpha, -beta, -gamma)
        + complex< double >(8.423, 0)       * Wigner::wigner_D_l2normalized(2,  0,  2, -alpha, -beta, -gamma)
        + complex< double >(1.8312, -9.8372)* Wigner::wigner_D_l2normalized(2,  0,  0, -alpha, -beta, -gamma);
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
                            grid(i, j, k) += complex< double >(cnt, cnt+1) * Wigner::wigner_D_l2normalized(l, M, Mp, -(2*M_PI*j)/(2*B), -(M_PI * (2.0 * k + 1.0)) / (4.0 * B), -(2*M_PI*i)/(2*B));
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
    DSOFTFourierCoefficients coef(bandwidth);
    
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
    DSOFTFourierCoefficients coef(bandwidth);
    DSOFTFourierCoefficients rec_coef(bandwidth);
    
    uniform_real_distribution< double > ctx;
    ctx.min = -1;
    ctx.max = +1;
    
    rand(coef, ctx);
    
    stopwatch sw = stopwatch::tic();
    FourierTransforms::IDSOFT(coef, sample);
    double time2 = sw.toc();
    
    // perform forward SOFT transform
    sw = stopwatch::tic();
    FourierTransforms::DSOFT(sample, rec_coef);
    double time = sw.toc();
    
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

void test_sofft_for_bandwidth(int B)
{
    grid3D< complex< double > > grid(2 * B, 0);
    DSOFTFourierCoefficients fc(B);
    
    // filling grid with custom wigners given in f
    for (int i = 0; i < 2 * B; ++i)
    {
        for (int j = 0; j < 2 * B; ++j)
        {
            for (int k = 0; k < 2 * B; ++k)
            {
                grid(i, j, k) = f((2*M_PI*j)/(2*B), (M_PI * (2.0 * k + 1.0)) / (4.0 * B), (2*M_PI*i)/(2*B));
            }
        }
    }
    
    FourierTransforms::DSOFT(grid, fc);
    
    printf("Spatial to spectral (DSOFT) for Bandwidth %d.\n\n", B);
    
    int cnt = 0;
    
    for (int i = 0; i < B; ++i)
    {
        if (cnt == 1000)
            break;
        
        for (int j = -i; j <= i; ++j)
        {
            if (cnt == 1000)
                break;
            
            for (int k = -i; k <= i; ++k)
            {
                if (cnt == 1000)
                    break;
                
                // Here the coefficients are printed out on the console
                printf("l=%4d, M=%4d, M'=%4d:    %s%.4f%s%.4f\n", i,  j, k, (fc(i,j,k).re < 0 ? "-" : " "), std::abs(fc(i,j,k).re), (fc(i,j,k).im < 0 ? " - " : " + "), std::abs(fc(i,j,k).im));
                
                cnt++;
            }
        }
    }
}

int main(int argc, const char ** argv)
{
    if (argc < 2)
    {
        printf("usage: ./soft_test <Bandwidth>\n");
        return 1;
    }
    
//    matrix< double > A(5, 5);
//    
//    A << 0.8147 << 0.0975 << 0.1576 << 0.1419 << 0.6557
//      << 0.9058 << 0.2785 << 0.9706 << 0.4218 << 0.0357
//      << 0.1270 << 0.5469 << 0.9572 << 0.9157 << 0.8491
//      << 0.9134 << 0.9575 << 0.4854 << 0.7922 << 0.9340
//      << 0.6324 << 0.9649 << 0.8003 << 0.9595 << 0.6787;
//    
//    vector< complex< double > > B(5, vec_type::COLUMN);
//    
//    B[0] = complex< double >(0.7577, 0.7060);
//    B[1] = complex< double >(0.7431, 0.0318);
//    B[2] = complex< double >(0.3922, 0.2769);
//    B[3] = complex< double >(0.6555, 0.0462);
//    B[4] = complex< double >(0.1712, 0.0971);
//    
//    std::cout << "A = " << A << std::endl;
//    std::cout << "B = " << B << std::endl;
//    std::cout << "-B = " << B * -1 << std::endl;
//    std::cout << "A * B = " << A * B << std::endl;
//    
//    B *= -1;
//    
//    std::cout << "B *= -1 = " << B << std::endl;
    
    matrix< complex< int > > A(5, 5);
    
    int B = atoi(*(argv + 1));
    for_back(B, false);
    
//    matrix< double > wig = DWT::wigner_d_matrix(5, 1, 2);
//    
//    std::cout << "wig = " << wig << std::endl;
    
    // create grid
//    test_sofft_for_bandwidth(B);
    
//    normal_distribution< double > nd;
//    nd.engine = random_engine::MERSENNE_TWISTER;
//    nd.mean = 0;
//    nd.standard_deviation = 1;
//    
//    uniform_int_distribution< int > uid;
//    uid.engine = random_engine::MERSENNE_TWISTER;
//    uid.min = 1;
//    uid.max = 4;
    
//    uniform_real_distribution< double > urd;
//    urd.engine = random_engine::MERSENNE_TWISTER64;
//    urd.min = -2;
//    urd.max = 2.342;
//    
//    uniform_int_distribution< int > uid;
//    uid.engine = random_engine::MERSENNE_TWISTER64;
//    uid.min = -1;
//    uid.max = +4;
//    
//    DSOFTFourierCoefficients fc(3);
//    
//    rand(fc, urd);
    
//    int l = 3, m = 0;
//    double x = constants< double >::pi / 6;
//    double y = constants< double >::pi / 3;
//    
//    stopwatch sw = stopwatch::tic();
//    complex< double > result = SphericalHarmonics::Ylm(l, m, x, y);
//    double time = sw.toc_micros();
//    
//    std::cout << std::fixed << result.re << (result.im < 0 ? " - " : " + ") << std::abs(result.im) << " in " << time << "micros" << std::endl;
    
    return 0;
}

#endif
