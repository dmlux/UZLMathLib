//
//  blas_wrapper.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 25.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_blas_wrapper_hpp
#define UZLMathLib_blas_wrapper_hpp

UZLMATH_BEGIN

/*- BLAS ENUM CUSTOMIZATIONS -*/
typedef enum UZLBLAS_ORDER     {UZLblasRowMajor, UZLblasColMajor}                                   UZLBLAS_ORDER;
typedef enum UZLBLAS_TRANSPOSE {UZLblasNoTrans, UZLblasTrans, UZLblasConjTrans, UZLblasConjNoTrans} UZLBLAS_TRANSPOSE;
typedef enum UZLBLAS_UPLO      {UZLblasUpper, UZLblasLower}                                         UZLBLAS_UPLO;
typedef enum UZLBLAS_DIAG      {UZLblasNonUnit, UZLblasUnit}                                        UZLBLAS_DIAG;
typedef enum UZLBLAS_SIDE      {UZLblasLeft, UZLblasRight}                                          UZLBLAS_SIDE;

/*- LAPACK CUSTOMIZATIONS -*/
#define uzlmath_complex_double double _Complex

extern "C"
{
    /*- CBLAS FUNCTIONS -*/

    void uzl_blas_cgemm(enum UZLBLAS_TRANSPOSE TransA, enum UZLBLAS_TRANSPOSE TransB, int M, int N, int K, float *alpha, float *A, int lda, float *B, int ldb, float *beta, float *C, int ldc);
    void uzl_blas_zgemm(enum UZLBLAS_TRANSPOSE TransA, enum UZLBLAS_TRANSPOSE TransB, int M, int N, int K, double *alpha, double *A, int lda, double *B, int ldb, double *beta, double *C, int ldc);
    void uzl_blas_sgemm(enum UZLBLAS_TRANSPOSE TransA, enum UZLBLAS_TRANSPOSE TransB, int M, int N, int K, float alpha, float *A, int lda, float *B, int ldb, float beta, float *C, int ldc);
    void uzl_blas_dgemm(enum UZLBLAS_TRANSPOSE TransA, enum UZLBLAS_TRANSPOSE TransB, int M, int N, int K, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc);
    
    /*- LAPACK FUNCTIONS -*/
    
    void uzl_lapack_dgetrf( int* m, int* n, double* a, int* lda, int* ipiv, int *info );
    void uzl_lapack_zgetrf( int* m, int* n, uzlmath_complex_double* a, int* lda, int* ipiv, int *info );
    
    /*- FFTW FUNCTIONS -*/
    
    void uzl_fftw_DFT2(size_t cols, size_t rows, double* arr);
    void uzl_fftw_layer_wise_DFT2_grid3D(int cols, int rows, int lays, double* arr);
    
    void uzl_fftw_IDFT2(size_t cols, size_t rows, double* arr);
    void uzl_fftw_layer_wise_IDFT2_grid3D(int cols, int rows, int lays, double* arr);
    
    void uzl_fftw_DFT(size_t size, double* arr);
    void uzl_fftw_IDFT(size_t size, double* arr);
}

UZLMATH_END

#endif
