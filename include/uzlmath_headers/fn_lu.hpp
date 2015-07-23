//
//  fn_matrix_lu.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 17.01.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_fn_matrix_lu_hpp
#define uzlmath_fn_matrix_lu_hpp

UZLMATH_BEGIN

/*! 
 * @brief   Calculates the LU-decomposition of a given matrix A.
 *
 * @since   0.0.1
 *
 * @author  Denis-Michael Lux <denis.lux@icloud.com>
 * @date    17.01.15
 */
template< typename eT >
inline
matrix< double > lu(const matrix< eT >& A)
{
    /*
     * INFO     = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
     *                has been completed, but the factor U is exactly
     *                singular, and division by zero will occur if it is used
     *                to solve a system of equations.
     */
    /*
    int INFO = (rows > cols) ? rows : cols;
    int LDA  = rows;
    
    int ipiv_dim = (rows < cols) ? rows : cols;
    int IPIV[ipiv_dim];
    
    int M = rows;
    int N = cols;
    
    matrix< double > lu(M, N);
    
    size_t i, cap = M * N;
    for (i = 0; i < cap; ++i)
    {
        lu.mem[i] = (double)mem[i];
    }
    
    LAPACK_dgetrf(&M, &N, lu.mem, &LDA, IPIV, &INFO);
    
    if ( INFO > 0 )
    {
        printf("** uzlmath error: LU-Decomposition argument error. Argument number %i was illegal. **\n", INFO);
        exit(EXIT_FAILURE);
    }
    
    return lu;
     */
    return matrix< double >();
}

template< typename eT >
inline
matrix< double > lu(const matrix< eT >& A, matrix< int >& P)
{
    /*
     * INFO     = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
     *                has been completed, but the factor U is exactly
     *                singular, and division by zero will occur if it is used
     *                to solve a system of equations.
     */
    /*
     int INFO = (rows > cols) ? rows : cols;
     int LDA  = rows;
     
     int ipiv_dim = (rows < cols) ? rows : cols;
     int IPIV[ipiv_dim];
     
     int M = rows;
     int N = cols;
     
     matrix< double > lu(M, N);
     size_t i, j, cap = M * N;
     for (i = 0; i < cap; ++i)
     {
     lu.mem[i] = (double)mem[i];
     }
     
     LAPACK_dgetrf(&M, &N, lu.mem, &LDA, IPIV, &INFO);
     
     if ( INFO > 0 )
     {
     printf("** uzlmath error: LU-Decomposition argument error. Argument number %i was illegal. **\n", INFO);
     exit(EXIT_FAILURE);
     }
     
     P = matrix< double >(M, M);
     P.eye();
     
     for (i = 0; i < ipiv_dim; ++i)
     {
     for (j = 0; j < P.n_cols(); ++j)
     {
     eT tmp = P(i, j);
     P(i, j) = P(IPIV[i] - 1, j);
     P(IPIV[i] - 1, j) = tmp;
     }
     }
     
     return lu;
     */
    return matrix< double >();
}

template< typename eT >
inline
void lu(const matrix< eT >& A, matrix< float >& LU, matrix< int >& P)
{
    /*
     * INFO     = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
     *                has been completed, but the factor U is exactly
     *                singular, and division by zero will occur if it is used
     *                to solve a system of equations.
     */
    /*
     int INFO = (rows > cols) ? rows : cols;
     int LDA  = rows;
     
     int ipiv_dim = (rows < cols) ? rows : cols;
     int IPIV[ipiv_dim];
     
     int M = rows;
     int N = cols;
     
     LU = matrix< float >(M, N);
     size_t i, j, cap = M * N;
     for (i = 0; i < cap; ++i)
     {
     LU.mem[i] = (float)mem[i];
     }
     
     LAPACK_sgetrf(&M, &N, LU.mem, &LDA, IPIV, &INFO);
     
     if ( INFO > 0 )
     {
     printf("** uzlmath error: LU-Decomposition argument error. Argument number %i was illegal. **\n", INFO);
     exit(EXIT_FAILURE);
     }
     
     P = matrix< float >(M, M);
     P.eye();
     
     for (i = 0; i < ipiv_dim; ++i)
     {
     for (j = 0; j < P.n_cols(); ++j)
     {
     eT tmp = P(i, j);
     P(i, j) = P(IPIV[i] - 1, j);
     P(IPIV[i] - 1, j) = tmp;
     }
     }
     */
}

template< typename eT >
inline
void lu(const matrix< eT >& A, matrix< double >& LU, matrix< int >& P)
{
    
}

template< typename eT >
inline
void lu(const matrix< eT >& A, matrix< float >& L , matrix< float >& U, matrix< double >& P)
{
    /*
     * INFO     = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
     *                has been completed, but the factor U is exactly
     *                singular, and division by zero will occur if it is used
     *                to solve a system of equations.
     */
    /*
     int INFO = (rows > cols) ? rows : cols;
     int LDA  = rows;
     
     int ipiv_dim = (rows < cols) ? rows : cols;
     int IPIV[ipiv_dim];
     
     int M = rows;
     int N = cols;
     
     double *lu = new double[M * N];
     
     size_t i, j, cap = M * N;
     for (i = 0; i < cap; ++i)
     {
     lu[i] = (double)mem[i];
     }
     
     LAPACK_dgetrf(&M, &N, lu, &LDA, IPIV, &INFO);
     
     if ( INFO > 0 )
     {
     printf("** uzlmath error: LU-Decomposition argument error. Argument number %i was illegal. **\n", INFO);
     exit(EXIT_FAILURE);
     }
     
     P = matrix< double >(M, M);
     P.eye();
     
     for (i = 0; i < ipiv_dim; ++i)
     {
     for (j = 0; j < P.n_cols(); ++j)
     {
     eT tmp = P(i, j);
     P(i, j) = P(IPIV[i] - 1, j);
     P(IPIV[i] - 1, j) = tmp;
     }
     }
     
     if (cols < rows)
     {
     L = matrix< double >(rows, cols);
     U = matrix< double >(cols, cols);
     }
     else
     {
     L = matrix< double >(rows, rows);
     U = matrix< double >(rows, cols);
     }
     
     for (i = 0; i < rows; ++i)
     {
     for (j = 0; j < cols; ++j)
     {
     if (j == i)
     {
     U(i, j) = lu[j * rows + i];
     L(i, j) = 1.0;
     }
     else if (j > i)
     {
     U(i, j) = lu[j * rows + i];
     }
     else
     {
     L(i, j) = lu[j * rows + i];
     }
     }
     }
     
     delete [] lu;
     */
}

template< typename eT >
inline
void lu(const matrix< eT >& A, matrix< double >& L , matrix< double >& U, matrix< double >& P)
{
    
}

UZLMATH_END

#endif
