//
//  fn_inv.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 17.01.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_fn_inv_hpp
#define UZLMathLib_fn_inv_hpp

UZLMATH_BEGIN

template< typename T >
inline
matrix< double > inv(const matrix< T >& A)
{
//    if ( A.n_rows() != A.n_cols() )
//    {
//        std::cout << "** uzlmath error: Try to invert a non-square matrix. **" << std::endl;
//        exit(EXIT_FAILURE);
//    }
//    
//    /*
//     * INFO     = 0:  successful exit
//     *          < 0:  if INFO = -i, the i-th argument had an illegal value
//     *          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
//     *                has been completed, but the factor U is exactly
//     *                singular, and division by zero will occur if it is used
//     *                to solve a system of equations.
//     */
//    int INFO = (A.n_rows() > A.n_cols()) ? A.n_rows() : A.n_cols();
//    int LDA  = A.n_rows();
//    
//    int ipiv_dim = (A.n_rows() < A.n_cols()) ? A.n_rows() : A.n_cols();
//    int IPIV[ipiv_dim];
//    
//    int M = A.n_rows();
//    int N = A.n_cols();
//    
//    matrix< double > inverse(M, N);
//    double* tmp_mem = inverse.memptr();
//    
//    if(is_double<eT>::value == true)
//    {
//        size_t cap = M * N;
//        memcpy(tmp_mem, A.memptr(), cap * sizeof(double));
//    }
//    else
//    {
//        const T* tmp_A = A.memptr();
//        size_t i, cap = M * N;
//        for(i = 0; i < cap; ++i)
//        {
//            tmp_mem[i] = (double)tmp_A[i];
//        }
//    }
//    
//    LAPACK_dgetrf(&M, &N, tmp_mem, &LDA, IPIV, &INFO);
//    
//    if ( INFO > 0 )
//    {
//        std::cout << "** uzlmath warning: Try to compute the inverse of singular matrix. **" << std::endl;
//        
//        size_t cap = M * N;
//        std::fill(tmp_mem, tmp_mem + cap, std::numeric_limits<double>::infinity());
//        
//        return inverse;
//    }
//    
//    int LWORK    = M * N;
//    double *WORK = new double[LWORK];
//    
//    LAPACK_dgetri(&M, tmp_mem, &M, IPIV, WORK, &LWORK, &INFO);
//    
//    delete [] WORK;
//    return inverse;
    return matrix< double >();
}

UZLMATH_END

#endif
