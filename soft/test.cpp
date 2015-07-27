
#include <iostream>
#include <uzlmath>

using namespace uzlmath;

int main(int argc, const char** argv)
{
    
    // creating matrices
    int n = 7, m = 6;
    matrix< double > A(n, m), C(m, n);
    
    // creating SOFT fourier coefficients for bandwidth B
    int B = 4;
    SOFTFourierCoefficients fc(B);
    SOFTFourierCoefficients fc_rec(B);
    
    // filling stuff with random values
    rand(fc, -1, 1);    // randoms between -1 and 1
    rand(A, -1, 1);     // randoms between -1 and 1
    rand(C, 1, 2);      // randoms between 1 and 2
    
    // printing matrices
    std::cout << "A = " << A << std::endl << "B = " << C << std::endl;
    std::cout << "C = A * B = " << A * C << std::endl;
    
    // creating grid for SOFT sample
    grid3D< complex< double > > sample(2 * B);
    
    // do inverse SOFT and measure time. After the inverse SOFT
    // the grid was filled with values that can be used for the
    // forward SOFT
    stopwatch sw = stopwatch::tic();
    FourierTransforms::ISOFT(fc, sample);
    double inv_time = sw.toc();
    
    // do forward SOFT and measure time
    sw = stopwatch::tic();
    FourierTransforms::SOFT(sample, fc_rec);
    double for_time = sw.toc();
    
    // print runtimes
    std::cout << "Inverse SOFT runtime: " << inv_time << "s" << std::endl;
    std::cout << "Forward SOFT runtime: " << for_time << "s" << std::endl << std::endl;
    
    // check if values are correct
    bool equal = true;
    double epsilon = 1e-13;
    for (int l = 0; l < B; ++l)
    {
        if ( !equal )
            break;
        
        for (int M = -l; M <= l; ++M)
        {
            if ( !equal )
                break;
            
            for (int Mp = -l; Mp <= l; ++Mp)
            {
                if (std::abs(fc(l, M, Mp).re - fc_rec(l, M, Mp).re) > epsilon && std::abs(fc(l, M, Mp).im - fc_rec(l, M, Mp).im) > epsilon)
                {
                    equal = false;
                    break;
                }
            }
        }
    }
    
    // print out other data
    std::cout << "Reconstructed coefficients and original equal? " << (equal ? "Yes" : "No") << std::endl << std::endl;
    
    if (B <= 4)
    {
        std::cout << "sample = " << sample << std::endl;
//        std::cout << "Random Fourier coefficients = " << fc << std::endl;
//        std::cout << "Reconstructed Fourier coefficients = " << fc_rec << std::endl;
    }
    
    return 0;
}