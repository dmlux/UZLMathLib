
/** CBLAS/LAPACK FROM OpenBLAS  **/
/** for linear Algebra puproses **/
#include <cblas.h>
#include <lapacke.h>

/** FFTW **/
/** for Fourier Analysis purposes **/
#include <fftw3.h>

// including function wrapper
#include <uzlmath>

#if _OPENMP
    #include <omp.h>
#endif

/*!
 * @brief           The namespace containg all includes of datatypes, functions and 
 *                  other namespaces for mathmatical computation from the **UZLMath**
 *                  library
 * @details         All datatypes and functions that are created for the library
 *                  'UZLMath' are encapsulated in this namespace.
 *
 * @since           0.0.1
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            14.05.2015
 */
namespace uzlmath
{
    /*- Wrapper for needed library functions -*/
    extern "C"
    {
        /*!
         * @brief           A collection of wrapper functions for external library functions.
         * @details         The functions are external functions from libraries. To provide
         *                  an independet usage of the **UZLMath** library the functions are
         *                  wrapped in own functions and then compiled into the **UZLMathLib**
         *                  library binary.
         * @defgroup        libwrapper Library wrapper functions
         * @{
         */

        void uzl_blas_cgemm(enum UZLBLAS_TRANSPOSE TransA, enum UZLBLAS_TRANSPOSE TransB, int M, int N, int K, float *alpha, float *A, int lda, float *B, int ldb, float *beta, float *C, int ldc)
        {
            // translate transposition
            CBLAS_TRANSPOSE ta, tb;
            if(TransA == UZLblasNoTrans)
            {
                ta = CblasNoTrans;
            }
            else if(TransA == UZLblasTrans)
            {
                ta = CblasTrans;
            }
            else if(TransA == UZLblasConjTrans)
            {
                ta = CblasConjTrans;
            }
            else
            {
                ta = CblasConjNoTrans;
            }
            
            // translate transposition
            if(TransB == UZLblasNoTrans)
            {
                tb = CblasNoTrans;
            }
            else if(TransB == UZLblasTrans)
            {
                tb = CblasTrans;
            }
            else if(TransB == UZLblasConjTrans)
            {
                tb = CblasConjTrans;
            }
            else
            {
                tb = CblasConjNoTrans;
            }
            
            cblas_cgemm(CblasColMajor, ta, tb, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
        }
        
        void uzl_blas_zgemm(enum UZLBLAS_TRANSPOSE TransA, enum UZLBLAS_TRANSPOSE TransB, int M, int N, int K, double *alpha, double *A, int lda, double *B, int ldb, double *beta, double *C, int ldc)
        {
            // translate transposition
            CBLAS_TRANSPOSE ta, tb;
            if(TransA == UZLblasNoTrans)
            {
                ta = CblasNoTrans;
            }
            else if(TransA == UZLblasTrans)
            {
                ta = CblasTrans;
            }
            else if(TransA == UZLblasConjTrans)
            {
                ta = CblasConjTrans;
            }
            else
            {
                ta = CblasConjNoTrans;
            }
            
            // translate transposition
            if(TransB == UZLblasNoTrans)
            {
                tb = CblasNoTrans;
            }
            else if(TransB == UZLblasTrans)
            {
                tb = CblasTrans;
            }
            else if(TransB == UZLblasConjTrans)
            {
                tb = CblasConjTrans;
            }
            else
            {
                tb = CblasConjNoTrans;
            }
            
            cblas_zgemm(CblasColMajor, ta, tb, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
        }
        
        void uzl_blas_sgemm(enum UZLBLAS_TRANSPOSE TransA, enum UZLBLAS_TRANSPOSE TransB, int M, int N, int K, float alpha, float *A, int lda, float *B, int ldb, float beta, float *C, int ldc)
        {
            // translate transposition
            CBLAS_TRANSPOSE ta, tb;
            if(TransA == UZLblasNoTrans)
            {
                ta = CblasNoTrans;
            }
            else if(TransA == UZLblasTrans)
            {
                ta = CblasTrans;
            }
            else if(TransA == UZLblasConjTrans)
            {
                ta = CblasConjTrans;
            }
            else
            {
                ta = CblasConjNoTrans;
            }
            
            // translate transposition
            if(TransB == UZLblasNoTrans)
            {
                tb = CblasNoTrans;
            }
            else if(TransB == UZLblasTrans)
            {
                tb = CblasTrans;
            }
            else if(TransB == UZLblasConjTrans)
            {
                tb = CblasConjTrans;
            }
            else
            {
                tb = CblasConjNoTrans;
            }
            
            cblas_sgemm(CblasColMajor, ta, tb, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
        }
        
        void uzl_blas_dgemm(enum UZLBLAS_TRANSPOSE TransA, enum UZLBLAS_TRANSPOSE TransB, int M, int N, int K, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc)
        {
            // translate transposition
            CBLAS_TRANSPOSE ta, tb;
            if(TransA == UZLblasNoTrans)
            {
                ta = CblasNoTrans;
            }
            else if(TransA == UZLblasTrans)
            {
                ta = CblasTrans;
            }
            else if(TransA == UZLblasConjTrans)
            {
                ta = CblasConjTrans;
            }
            else
            {
                ta = CblasConjNoTrans;
            }
            
            // translate transposition
            if(TransB == UZLblasNoTrans)
            {
                tb = CblasNoTrans;
            }
            else if(TransB == UZLblasTrans)
            {
                tb = CblasTrans;
            }
            else if(TransB == UZLblasConjTrans)
            {
                tb = CblasConjTrans;
            }
            else
            {
                tb = CblasConjNoTrans;
            }
            
            cblas_dgemm(CblasColMajor, ta, tb, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
        }
        
        
        
        
        void uzl_lapack_dgetrf( int* m, int* n, double* a, int* lda, int* ipiv, int *info )
        {
            LAPACK_dgetrf( m, n, a, lda, ipiv, info );
        }
        
        void uzl_lapack_zgetrf( int* m, int* n, uzlmath_complex_double* a, int* lda, int* ipiv, int *info )
        {
            LAPACK_zgetrf( m, n, a, lda, ipiv, info );
        }
        
        
        
        void uzl_fftw_layer_wise_fft2_cube(int cols, int rows, int lays, double* arr)
        {
            #if _OPENMP
            
            // define indices
            unsigned int i;
            
            // storage plans
            fftw_plan plans[lays];
            
            // create plans
            // Creating plans and destroying them is not really thread safe
            // therefore the plans should be created and destoryed sequentially
            for (i = 0; i < lays; ++i)
            {
                // get correct layer
                fftw_complex* layer = (fftw_complex*)arr + i * rows * cols;
                
                // create plan to execute an FFT
                plans[i] = fftw_plan_dft_2d(cols, rows, layer, layer, FFTW_FORWARD, FFTW_ESTIMATE);
            }
            
            // execute plans in parallel
            #pragma omp parallel for private(i) shared(lays, plans) schedule(dynamic)
            for (i = 0; i < lays; ++i)
            {
                // execute FFT2 plan
                fftw_execute(plans[i]);
            }
            
            for (i = 0; i < lays; ++i)
            {
                fftw_destroy_plan(plans[i]);
            }
            
            // remove all additional allocated objects needed for the FFT
            fftw_cleanup();
            
            #else
            
            // define parameters for many dft
            int rank    = 2;            // Dimension of FFT's
            int n[]     = {cols, rows}; // Dimension of matrix for each FFT
            int howmany = lays;         // How many FFT's should be executed?
            int idist   = n[0] * n[1];  // Length of memory space for each matrix
            int odist   = idist;        // ... for output same, but is in-place...
            int istride = 1;            // 1 because matrices are contigous in memory
            int ostride = 1;            // ... for output same, but is still in-place...
            int* inembed= n;
            int* onembed= n;
            
            // create plan to execute an layerwise FFT2
            fftw_plan lay_wise_fft2 = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*)arr, inembed, istride, idist, (fftw_complex*)arr, onembed, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);
            
            // execute layer-wise FFT2
            fftw_execute(lay_wise_fft2);
            
            // free allocated memory
            fftw_destroy_plan(lay_wise_fft2);
            fftw_cleanup();
            
            #endif
        }
        
        void uzl_fftw_layer_wise_ifft2_cube(int cols, int rows, int lays, double*  arr)
        {
            #if _OPENMP
            
            // define indices
            unsigned int i;
            
            // storage plans
            fftw_plan plans[lays];
            
            // create plans
            // Creating plans and destroying them is not really thread safe
            // therefore the plans should be created and destoryed sequentially
            for (i = 0; i < lays; ++i)
            {
                // get correct layer
                fftw_complex* layer = (fftw_complex*)arr + i * rows * cols;
                
                // create plan to execute an FFT
                plans[i] = fftw_plan_dft_2d(cols, rows, layer, layer, FFTW_BACKWARD, FFTW_ESTIMATE);
            }
            
            // execute plans in parallel
            #pragma omp parallel for private(i) shared(lays, plans) schedule(dynamic)
            for (i = 0; i < lays; ++i)
            {
                // execute FFT2 plan
                fftw_execute(plans[i]);
            }
            
            for (i = 0; i < lays; ++i)
            {
                fftw_destroy_plan(plans[i]);
            }
            
            // remove all additional allocated objects needed for the FFT
            fftw_cleanup();
            
            #else
            
            // define parameters for many dft
            int rank    = 2;            // Dimension for IFFT's
            int n[]     = {cols, rows}; // Dimension of matrix for each IFFT
            int howmany = lays;         // How many IFFT's should be executed?
            int idist   = n[0] * n[1];  // Length of memory space for each matrix
            int odist   = idist;        // ... for output same, but is in-place...
            int istride = 1;            // 1 because matrices are contigous in memory
            int ostride = 1;            // ... for output same, but is still in-place...
            int* inembed= n;
            int* onembed= n;
            
            // create plan to execute an layerwise IFFT2
            fftw_plan lay_wise_ifft2 = fftw_plan_many_dft(rank, n, howmany, (fftw_complex*)arr, inembed, istride, idist, (fftw_complex*)arr, onembed, ostride, odist, FFTW_BACKWARD, FFTW_ESTIMATE);
            
            // execute layer-wise IFFT2
            fftw_execute(lay_wise_ifft2);
            
            // free allocated memory
            fftw_destroy_plan(lay_wise_ifft2);
            fftw_cleanup();
            
            #endif
        }
        
        void uzl_fftw_fft2(size_t cols, size_t rows, double* arr)
        {
            // allocate storage for the FFT
            //fftw_complex  inFFT2[cols * rows];
            fftw_complex outFFT2[cols * rows];
            
            // create plan to execute an FFT
            fftw_plan fft2 = fftw_plan_dft_2d(cols, rows, (fftw_complex*)arr, outFFT2, FFTW_FORWARD, FFTW_ESTIMATE);
            
            // execute FFT2 plan
            fftw_execute(fft2);
            
            memcpy(arr, outFFT2, rows * cols * sizeof(fftw_complex));
            
            // free allocated memory
            fftw_destroy_plan(fft2);
            
            // remove all additional allocated objects needed for the FFT
            fftw_cleanup();
        }
        
        void uzl_fftw_ifft2(size_t cols, size_t rows, double* arr)
        {
            // allocate storage for the FFT
            fftw_complex outFFT2[cols * rows];
            
            // create plan to execute an FFT
            // Since this is a C code section the double
            // array can simply be casted into fftw_complex.
            fftw_plan fft2 = fftw_plan_dft_2d(cols, rows, (fftw_complex*)arr, outFFT2, FFTW_BACKWARD, FFTW_ESTIMATE);
            
            // execute FFT2 plan
            fftw_execute(fft2);
            
            memcpy(arr, outFFT2, rows * cols * sizeof(fftw_complex));
            
            // free allocated memory
            fftw_destroy_plan(fft2);
            
            // remove all additional allocated objects needed for the FFT
            fftw_cleanup();
        }
        
        void uzl_fftw_fft(size_t size, double* arr)
        {
            // create FFT plan
            fftw_plan p = fftw_plan_dft_1d(size, (fftw_complex*)arr, (fftw_complex*)arr, FFTW_FORWARD, FFTW_ESTIMATE);
            
            // execute plan
            fftw_execute(p);
            
            // destroy plan and clean up workspace
            fftw_destroy_plan(p);
            fftw_cleanup();
        }
        
        void uzl_fftw_ifft(size_t size, double* arr)
        {
            // create IFFT plan
            fftw_plan p = fftw_plan_dft_1d(size, (fftw_complex*)arr, (fftw_complex*)arr, FFTW_BACKWARD, FFTW_ESTIMATE);
            
            // execute plan
            fftw_execute(p);
            
            // destroy plan and clean up workspace
            fftw_destroy_plan(p);
            fftw_cleanup();
        }
        
        /*!
         * @}
         */
    }
}