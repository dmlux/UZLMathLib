//
//  cx_matrix_def.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 13.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_cx_matrix_def_hpp
#define uzlmath_cx_matrix_def_hpp

UZLMATH_BEGIN

/*!
 * @brief           The deconstructor to delete a complex matrix and its contents properly from
 *                  RAM.
 * @details         Deletes the complex matrix and frees its allocated dynamic mamemory.
 */
template< typename eT >
inline
matrix< complex< eT > >::~matrix()
{
    delete [] mem;
}

/*!
 * @brief           Default constructor to build a complex matrix with no rows and columns.
 * @details         Constructs a complex matrix that has no columns and rows and is declared
 *                  as not initialized.
 */
template< typename eT >
inline
matrix< complex< eT > >::matrix()
    : rows(0)
    , cols(0)
    , r_inj(0)
    , c_inj(0)
    , mem(nullptr)
{}

/*!
 * @brief           Constructs a complex matrix with given parameters
 * @details         Constructs the complex matrix by allocating enough memory to store \f$M\times N\f$
 *                  values. The content of the complex matrix is initialized with \f$0 + 0i\f$
 *
 * @param[in]       m Number of rows in the constructed complex matrix
 * @param[in]       n Number of columns in the constructed complex matrix
 */
template< typename eT >
inline
matrix< complex< eT > >::matrix(const size_t& m, const size_t& n)
    : rows(m)
    , cols(n)
    , r_inj(0)
    , c_inj(0)
{
    size_t cap  = m * n;
    mem         = new complex< eT >[cap];
}

/*!
 * @brief           Constructs a square complex matrix with given dimension \f$N\times N\f$
 * @details         Constructs the complex matrix by allocating enough memory to store \f$N\times N\f$
 *                  values. The content of the complex matrix is initialized with \f$0+0i\f$.
 *
 * @param[in]       mn Number of rows and columns in the constructed complex matrix
 */
template< typename eT >
inline
matrix< complex< eT > >::matrix(const size_t& mn)
    : rows(mn)
    , cols(mn)
    , r_inj(0)
    , c_inj(0)
{
    size_t cap  = mn * mn;
    mem         = new complex< eT >[cap];
}

/*!
 * @brief           Constructs a complex matrix with given dimensions and given intial element value
 * @details         Constructs the complex matrix by allocating enough memory to store \f$M\times N\f$
 *                  values. Each value in the complex matrix will be initialized with given initial value.
 *
 * @param[in]       m Number of rows in the constructed complex matrix
 * @param[in]       n Number of columns in the constructed complex matrix
 * @param[in]       initial Initial value for each complex matrix element
 */
template< typename eT >
inline
matrix< complex< eT > >::matrix(const size_t& m, const size_t& n, const complex< eT >& initial)
    : rows(m)
    , cols(n)
    , r_inj(0)
    , c_inj(0)
{
    size_t cap  = m * n;
    mem         = new complex< eT >[cap];
    
    if (cap > 0)
    {
        for (int i = 0; i < cap; ++i)
        {
            mem[i] = initial;
        }
    }
}

/*!
 * @brief           A copy constructor for copying a given complex matrix
 * @details         Copies the contents of the given complex matrix to build a new complex matrix
 *                  with the same state.
 *
 * @param[in]       A The complex matrix that is supposed to be copied.
 */
template< typename eT >
inline
matrix< complex< eT > >::matrix(const matrix< complex< eT > >& A)
    : rows(A.rows)
    , cols(A.cols)
    , r_inj(A.r_inj)
    , c_inj(A.c_inj)
{
    size_t cap  = rows * cols;
    mem         = new complex< eT >[cap];
    
    if (cap > 0)
    {
        memcpy(mem, A.mem, cap * sizeof(complex< eT >));
    }
}

/*!
 * @brief           A copy constructor for copying a given matrix
 * @details         Copies the contents of the given matrix to build a new complex matrix
 *                  with the same state.
 *
 * @param[in]       A The matrix that is supposed to be copied.
 */
template< typename eT >
inline
matrix< complex< eT > >::matrix(const matrix< eT >& A)
    : rows(A.rows)
    , cols(A.cols)
    , r_inj(A.r_inj)
    , c_inj(A.c_inj)
{
    size_t cap  = rows * cols;
    mem         = new complex< eT >[cap];
    
    for (int i = 0; i < cap; ++i)
    {
        mem[i] = complex< eT >(A.mem[i], 0);
    }
}

/*!
 * @brief           A move constructor to copy an r-value complex matrix into a new complex matrix
 * @details         Takes the memory of the r-value complex matrix \f$A\f$ and copies all other
 *                  instance variables. The resulting comple matrix is in the same state than the 
 *                  original complex matrix. A move constructor has better performance by assigning 
 *                  complex matrices that were calculated in one expression.
 *
 * @param[in,out]   A The r-value complex matrix \f$A\f$ which content should be moved to a
 *                  new complex matrix.
 */
template< typename eT >
inline
matrix< complex< eT > >::matrix(matrix< complex< eT > >&& A)
    : rows(A.rows)
    , cols(A.cols)
    , r_inj(A.r_inj)
    , c_inj(A.c_inj)
{
    complex< eT >* tmp = mem;
    mem                = A.mem;
    A.mem              = tmp;
}



/*!
 * @brief           The element-wise addition operator for complex matrices.
 * @details         Adds the complex matrix \f$A\f$ to the current complex matrix by adding
 *                  the complex matrices element-wise. The complex matrices should have the
 *                  same dimensions to be added.
 *
 * @param[in]       A The complex matrix \f$A\f$ on the right handside of the addition operator
 *
 * @return          A new complex matrix containing the result of the addition.
 */
template< typename eT >
inline
matrix< complex< eT > > matrix< complex< eT > >::operator+(const matrix< complex< eT > >& A)
{
    if ( rows != A.rows || cols != A.cols )
    {
        uzlmath_error("%s", "Dimension mismatch in matrix-matrix addition.");
    }
    
    matrix< complex< eT > > C(rows, cols);
    size_t i, j;
    
    for (i = 0; i < rows; ++i)
    {
        for (j = 0; j < cols; ++j)
        {
            C(i, j) = mem[j * rows + i] + A(i, j);
        }
    }
    
    return C;
}

/*!
 * @brief           The element-wise addition operator for a non-complex matrices.
 * @details         Adds the given matrix \f$A\f$ to the current complex matrix by
 *                  adding both matrices element-wise. The given matrices should
 *                  have the same dimensions to be added.
 *
 * @param[in]       A The matrix \f$A\f$ on the right handside of the addition operator
 *
 * @return          A new complex matrix containing the result of the addition.
 */
template< typename eT >
inline
matrix< complex< eT > > matrix< complex< eT > >::operator+(const matrix< eT >& A)
{
    if ( rows != A.rows || cols != A.cols )
    {
        uzlmath_error("%s", "Dimension mismatch in matrix-matrix addition.");
    }
    
    matrix< complex< eT > > C(rows, cols);
    size_t i, j;
    
    for (i = 0; i < rows; ++i)
    {
        for (j = 0; j < cols; ++j)
        {
            C(i, j) = mem[j * rows + i] + complex< eT >(A(i, j), 0);
        }
    }
    
    return C;
}

/*!
 * @brief           The element-wise minus operator for two complex matrices.
 * @details         Subtracts the complex matrix \f$A\f$ from the current complex matrix by subtracting the
 *                  complex matrix elements of A element-wise. The complex matrices should have the same
 *                  dimensions to be subtracted
 *
 * @param[in]       A The complex matrix \f$A\f$ on the right hand-side of the minus operator
 *
 * @return          A new complex matrix containing the result of the subtraction
 */
template< typename eT >
inline
matrix< complex< eT > > matrix< complex< eT > >::operator-(const matrix< complex< eT > >& A)
{
    if ( rows != A.rows || cols != A.cols )
    {
        uzlmath_error("%s", "Dimension mismatch in matrix-matrix subtraction.");
    }
    
    matrix< complex< eT > > C(rows, cols);
    size_t i, j;
    
    for (i = 0; i < rows; ++i)
    {
        for (j = 0; j < cols; ++j)
        {
            C(i, j) = mem[j * rows + i] - A(i, j);
        }
    }
    
    return C;
}

/*!
 * @brief           The element-wise minus operator for non-complex matrix.
 * @details         Subtracts the given matrix \f$A\f$ from the current complex matrix by subtracting
 *                  the matrix elements of A element-wise. The given matrices should have the same
 *                  dimensions to be subtracted
 *
 * @param[in]       A The matrix \f$A\f$ on the right hand-side of the minus operator
 *
 * @return          A new complex matrix containing the result of the subtraction
 */
template< typename eT >
inline
matrix< complex< eT > > matrix< complex< eT > >::operator-(const matrix< eT >& A)
{
    if ( rows != A.rows || cols != A.cols )
    {
        uzlmath_error("%s", "Dimension mismatch in matrix-matrix subtraction.");
    }
    
    matrix< complex< eT > > C(rows, cols);
    size_t i, j;
    
    for (i = 0; i < rows; ++i)
    {
        for (j = 0; j < cols; ++j)
        {
            C(i, j) = mem[j * rows + i] - complex< eT >(A(i, j), 0);
        }
    }
    
    return C;
}

/*!
 * @brief           The multiplication operator for complex matrices.
 * @details         Multiplies the current complex matrix by the given matrix \f$A\f$ by applying
 *                  \f{eqnarray*}{
 *                      (CA)_{ij} = \sum\limits_{k=1}^mA_{ik}B_{kj}
 *                  \f}
 *                  Where \f$C\f$ (current matrix) is an \f$n\times m\f$ complex matrix, \f$A\f$ is
 *                  an \f$m\times p\f$ complex matrix and each \f$i,j\f$ entry is given by multiplying
 *                  the entries \f$A_{ik}\f$ (across row \f$i\f$ of \f$C\f$) by the entries
 *                  \f$B_{kj}\f$ (down column \f$j\f$) of \f$A\f$, for \f$k=0,1,\dots,m\f$,
 *                  and summing the results over \f$k\f$
 *
 * @param[in]       A The complex matrix that is supposed to be mutliplied with the current complex matrix
 *
 * @return          A new complex matrix containing the result of the multiplication
 */
template< typename eT >
inline
matrix< complex< eT > > matrix< complex< eT > >::operator*(const matrix< complex< eT > >& A)
{
    if ( cols != A.rows )
    {
        uzlmath_error("%s", "Dimension mismatch in matrix-matrix multiplication.");
    }
    
    size_t n_rows = rows;
    size_t n_cols = A.cols;
    
    matrix< complex< eT > > result(n_rows, n_cols);
    
    if (is_int< eT >::value == true || is_short< eT >::value == true)
    {
        
        float* tmp_mem  = new float[2 * rows * cols];
        float* tmp_A    = new float[2 * A.rows * A.cols];
        
        size_t i, cap = rows * cols, cap_a = A.rows * A.cols;
        for (i = 0; i < cap; ++i)
        {
            tmp_mem[i * 2]      = static_cast< float >(mem[i].re);
            tmp_mem[i * 2 + 1]  = static_cast< float >(mem[i].im);
        }
        
        for (i = 0; i < cap_a; ++i)
        {
            tmp_A[i * 2]        = static_cast< float >(A.mem[i].re);
            tmp_A[i * 2 + 1]    = static_cast< float >(A.mem[i].im);
        }
        
        float* C        = new float[2 * n_rows * n_cols];
        float alpha[2]  = {1, 0};
        float beta[2]   = {0, 0};
        uzlblas_cgemm(UZLblasNoTrans, UZLblasNoTrans, n_rows, n_cols, cols, alpha, tmp_mem, rows, tmp_A, A.rows, beta, C, n_rows);
        
        size_t cap_c = n_rows * n_cols;
        
        for (i = 0; i < cap_c; ++i)
        {
            result.mem[i] = complex< eT >(C[i * 2], C[i * 2 +1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_A;
        delete [] C;
    }
    else if (is_double< eT >::value == true)
    {
        // Treat pointers as double pointers.
        double* A_mem_ptr = reinterpret_cast< double* >(mem);
        double* B_mem_ptr = reinterpret_cast< double* >(A.mem);
        double* C_mem_ptr = reinterpret_cast< double* >(result.mem);
        
        double alpha[2]   = {1, 0};
        double beta[2]    = {0, 0};
        uzlblas_zgemm(UZLblasNoTrans, UZLblasNoTrans, n_rows, n_cols, cols, alpha,
                      A_mem_ptr, rows, B_mem_ptr, A.rows, beta, C_mem_ptr, n_rows);
    }
    else if (is_float< eT >::value == true)
    {
        // Treat pointers as float pointers.
        float* A_mem_ptr = reinterpret_cast< float* >(mem);
        float* B_mem_ptr = reinterpret_cast< float* >(A.mem);
        float* C_mem_ptr = reinterpret_cast< float* >(result.mem);
        
        float alpha[2]   = {1, 0};
        float beta[2]    = {0, 0};
        uzlblas_cgemm(UZLblasNoTrans, UZLblasNoTrans, n_rows, n_cols, cols, alpha,
                      A_mem_ptr, rows, B_mem_ptr, A.rows, beta, C_mem_ptr, n_rows);
    }
    else
    {
        double* tmp_mem = new double[2 * rows * cols];
        double* tmp_A   = new double[2 * A.rows * A.cols];
        
        size_t i, cap = rows * cols, cap_a = A.rows * A.cols;
        for (i = 0; i < cap; ++i)
        {
            tmp_mem[i * 2]      = static_cast< double >(mem[i].re);
            tmp_mem[i * 2 + 1]  = static_cast< double >(mem[i].im);
        }
        
        for (i = 0; i < cap_a; ++i)
        {
            tmp_A[i * 2]        = static_cast< double >(A.mem[i].re);
            tmp_A[i * 2 + 1]    = static_cast< double >(A.mem[i].im);
        }
        
        double* C       = new double[2 * n_rows * n_cols];
        double alpha[2] = {1, 0};
        double beta[2]  = {0, 0};
        uzlblas_zgemm(UZLblasNoTrans, UZLblasNoTrans, n_rows, n_cols, cols, alpha, tmp_mem, rows, tmp_A, A.rows, beta, C, n_rows);
        
        size_t cap_c = n_rows * n_cols;
        
        for (i = 0; i < cap_c; ++i)
        {
            result.mem[i] = complex< eT >(C[i * 2], C[i * 2 + 1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_A;
        delete [] C;
    }
    
    return result;
}

/*!
 * @brief           The multiplication operator for a non-complex matrix.
 * @details         Multiplies the current complex matrix by the given matrix \f$A\f$ by applying
 *                  \f{eqnarray*}{
 *                      (CA)_{ij} = \sum\limits_{k=1}^mA_{ik}B_{kj}
 *                  \f}
 *                  Where \f$C\f$ (current matrix) is an \f$n\times m\f$ complex matrix, \f$A\f$ is
 *                  an \f$m\times p\f$ non-complex matrix and each \f$i,j\f$ entry is given by multiplying
 *                  the entries \f$A_{ik}\f$ (across row \f$i\f$ of \f$C\f$) by the entries
 *                  \f$B_{kj}\f$ (down column \f$j\f$) of \f$A\f$, for \f$k=0,1,\dots,m\f$,
 *                  and summing the results over \f$k\f$
 *
 * @param[in]       A The non-complex matrix that is supposed to be mutliplied with the current complex matrix
 *
 * @return          A new complex matrix containing the result of the multiplication
 */
template< typename eT >
inline
matrix< complex< eT > > matrix< complex< eT > >::operator*(const matrix< eT >& A)
{
    if ( cols != A.rows )
    {
        uzlmath_error("%s", "Dimension mismatch in matrix-matrix multiplication.");
    }
    
    size_t n_rows = rows;
    size_t n_cols = A.cols;
    
    matrix< complex< eT > > result(n_rows, n_cols);
    
    if (is_int< eT >::value == true || is_long< eT >::value == true || is_short< eT >::value == true)
    {
        
        float* tmp_mem  = new float[2 * rows * cols];
        float* tmp_A    = new float[2 * A.rows * A.cols];
        
        size_t i, cap = rows * cols, cap_a = A.rows * A.cols;
        for (i = 0; i < cap; ++i)
        {
            tmp_mem[i * 2]      = static_cast< float >(mem[i].re);
            tmp_mem[i * 2 + 1]  = static_cast< float >(mem[i].im);
        }
        
        for (i = 0; i < cap_a; ++i)
        {
            tmp_A[i * 2]        = static_cast< float >(A.mem[i]);
            tmp_A[i * 2 + 1]    = static_cast< float >(0);
        }
        
        float* C        = new float[2 * n_rows * n_cols];
        float alpha[2]  = {1, 0};
        float beta[2]   = {0, 0};
        uzlblas_cgemm(UZLblasNoTrans, UZLblasNoTrans, n_rows, n_cols, cols, alpha, tmp_mem, rows, tmp_A, A.rows, beta, C, n_rows);
        
        size_t cap_c = n_rows * n_cols;
        
        for (i = 0; i < cap_c; ++i)
        {
            result.mem[i] = complex< eT >(C[i * 2], C[i * 2 +1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_A;
        delete [] C;
    }
    else if (is_double< eT >::value == true)
    {
        double* tmp_mem = reinterpret_cast< double* >(mem);
        double* tmp_A   = new double[2 * A.rows * A.cols];
        
        size_t i, cap_a = A.rows * A.cols;
        for (i = 0; i < cap_a; ++i)
        {
            tmp_A[i * 2]        = static_cast< double >(A.mem[i]);
            tmp_A[i * 2 + 1]    = static_cast< double >(0);
        }
        
        double* C       = reinterpret_cast< double* >(result.mem);
        double alpha[2] = {1, 0};
        double beta[2]  = {0, 0};
        uzlblas_zgemm(UZLblasNoTrans, UZLblasNoTrans, n_rows, n_cols, cols, alpha, tmp_mem, rows, tmp_A, A.rows, beta, C, n_rows);
        
        delete [] tmp_A;
    }
    else if (is_float< eT >::value == true)
    {
        float* tmp_mem = reinterpret_cast< float* >(mem);
        float* tmp_A   = new float[2 * A.rows * A.cols];
        
        size_t i, cap_a = A.rows * A.cols;
        for (i = 0; i < cap_a; ++i)
        {
            tmp_A[i * 2]        = static_cast< float >(A.mem[i]);
            tmp_A[i * 2 + 1]    = static_cast< float >(0);
        }
        
        float* C       = reinterpret_cast< float* >(result.mem);
        float alpha[2] = {1, 0};
        float beta[2]  = {0, 0};
        uzlblas_cgemm(UZLblasNoTrans, UZLblasNoTrans, n_rows, n_cols, cols, alpha, tmp_mem, rows, tmp_A, A.rows, beta, C, n_rows);
        
        delete [] tmp_A;
    }
    else
    {
        double* tmp_mem = new double[2 * rows * cols];
        double* tmp_A   = new double[2 * A.rows * A.cols];
        
        size_t i, cap = rows * cols, cap_a = A.rows * A.cols;
        for (i = 0; i < cap; ++i)
        {
            tmp_mem[i * 2]      = static_cast< double >(mem[i].re);
            tmp_mem[i * 2 + 1]  = static_cast< double >(mem[i].im);
        }
        
        for (i = 0; i < cap_a; ++i)
        {
            tmp_A[i * 2]        = static_cast< double >(A.mem[i]);
            tmp_A[i * 2 + 1]    = static_cast< double >(0);
        }
        
        double* C       = new double[2 * n_rows * n_cols];
        double alpha[2] = {1, 0};
        double beta[2]  = {0, 0};
        uzlblas_zgemm(UZLblasNoTrans, UZLblasNoTrans, n_rows, n_cols, cols, alpha, tmp_mem, rows, tmp_A, A.rows, beta, C, n_rows);
        
        size_t cap_c = n_rows * n_cols;
        
        for (i = 0; i < cap_c; ++i)
        {
            result.mem[i] = complex< eT >(C[i * 2], C[i * 2 + 1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_A;
        delete [] C;
    }
    
    return result;
}

/*!
 * @brief           The assignment operator to copy matrices.
 * @details         Copies the given matrix \f$A\f$ into a new complex matrix
 *
 * @param[in]       A The matrix that is supposed to be copied
 *
 * @return          A new matrix that has the same state than \f$A\f$ and containing
 *                  the same values.
 */
template< typename eT >
inline
const matrix< complex< eT > >& matrix< complex< eT > >::operator=(const matrix< eT >& A)
{
    rows = A.rows;
    cols = A.cols;
    size_t cap = rows * cols;
    
    delete [] mem;
    mem = new complex< eT >[cap];
    
    if (cap > 0)
    {
        for (int i = 0; i < cap; ++i)
        {
            mem[i] = complex< eT >(A.mem[i], 0);
        }
    }
    
    return *this;
}

/*!
 * @brief           The assignment operator to copy a complex matrices.
 * @details         Copies the given complex matrix \f$A\f$ into a new complex matrix
 *
 * @param[in]       A The complex matrix that is supposed to be copied
 *
 * @return          A new complex matrix that has the same state than \f$A\f$ and containing
 *                  the same values.
 */
template< typename eT >
inline
const matrix< complex< eT > >& matrix< complex< eT > >::operator=(const matrix< complex< eT > >& A)
{
    if ( this == &A )
    {
        return *this;
    }
    
    rows = A.rows;
    cols = A.cols;
    size_t cap = rows * cols;
    
    delete [] mem;
    mem = new complex< eT >[cap];
    
    if (cap > 0)
    {
        memcpy(mem, A.mem, cap * sizeof(complex< eT >));
    }
        
    return *this;
}

/*!
 * @brief           The assignment operator to copy matrices by moving its contents.
 * @details         Moving the contents of a r-value matrix into a new matrix
 *
 * @param[in]       A The r-value matrix that is supposed to be copied
 *
 * @return          A new matrix that has the same state than \f$A\f$ and containing
 *                  the same values.
 */
template< typename eT >
inline
const matrix< complex< eT > >& matrix< complex< eT > >::operator=(matrix< complex< eT > >&& A)
{
    if ( this == &A )
    {
        return *this;
    }
    
    rows = A.rows;
    cols = A.cols;
    
    complex< eT >* tmp  = mem;
    mem                 = A.mem;
    A.mem               = tmp;
    
    return *this;
}



/*!
 * @brief           The _equals_ compare operator to compare two matrices element-wise.
 * @details         Comparing the current matrix with the matrix \f$A\f$ by checking each
 *                  element with the same indices in both matrices. Both matrices should
 *                  have the same dimensions.
 *
 * @param[in]       A The matrix that is supposed to be compared to the current matrix
 *
 * @return          True if both matrices are equal, false else.
 */
template< typename eT >
inline
bool matrix< complex< eT > >::operator==(const matrix< complex< eT > >& A)
{
    if ( rows != A.rows || cols != A.cols )
    {
        return false;
    }
    
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        if (mem[i] != A.mem[i])
        {
            return false;
        }
    }
    
    return true;
}

/*!
 * @brief           The _unequals_ compare operator to compare two matrices element-wise
 * @details         Comparing the current matrix with the matrix \f$A\f$ by checking each
 *                  element with the same indicies in both matrices. Both matrices should
 *                  have the same dimensions
 *
 * @param[in]       A The matrix that is supposed to be compared to the current matrix.
 *
 * @return          True if both matrices are unequal, false else.
 */
template< typename eT >
inline
bool matrix< complex< eT > >::operator!=(const matrix< complex< eT > >& A)
{
    return !(*this == A);
}



/*!
 * @brief           The addition assignment operator.
 * @details         Adds the matrix \f$A\f$ to the current matrix and stores the resulting
 *                  matrix in the current matrix. Both matrices should have the same size.
 *
 * @param[in]       A The matrix that is supposed to be added to the current matrix.
 *
 * @return          The reference to the current matrix that contains the result.
 */
template< typename eT >
inline
const matrix< complex< eT > >& matrix< complex< eT > >::operator+=(const matrix< eT >& A)
{
    if ( rows != A.rows || cols != A.cols )
    {
        uzlmath_error("%s", "Dimension mismatch in matrix-matrix addition.");
    }
    
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        mem[i] += complex< eT >(A.mem[i], 0);
    }
    
    return *this;
}

/*!
 * @brief           The addition assignment operator.
 * @details         Adds the matrix \f$A\f$ to the current matrix and stores the resulting
 *                  matrix in the current matrix. Both matrices should have the same size.
 *
 * @param[in]       A The matrix that is supposed to be added to the current matrix.
 *
 * @return          The reference to the current matrix that contains the result.
 */
template< typename eT >
inline
const matrix< complex< eT > >& matrix< complex< eT > >::operator+=(const matrix< complex< eT > >& A)
{
    if ( rows != A.rows || cols != A.cols )
    {
        uzlmath_error("%s", "Dimension mismatch in matrix-matrix addition.");
    }
    
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        mem[i] += A.mem[i];
    }
    
    return *this;
}

/*!
 * @brief           The subtraction assignment operator.
 * @details         Subtracts the matrix \f$A\f$ from the current matrix and stores the
 *                  result in the current matrix. Both matrices should have the same size
 *
 * @param[in]       A The matrix that is supposed to be subtracted from the current matrix
 *
 * @return          The reference to the current matrix that contains the result.
 */
template< typename eT >
inline
const matrix< complex< eT > >& matrix< complex< eT > >::operator-=(const matrix< eT >& A)
{
    if ( rows != A.rows || cols != A.cols )
    {
        uzlmath_error("%s", "Dimension mismatch in matrix-matrix subtraction.");
    }
    
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        mem[i] -= complex< eT >(A.mem[i], 0);
    }
    
    return *this;
}

/*!
 * @brief           The subtraction assignment operator.
 * @details         Subtracts the matrix \f$A\f$ from the current matrix and stores the
 *                  result in the current matrix. Both matrices should have the same size
 *
 * @param[in]       A The matrix that is supposed to be subtracted from the current matrix
 *
 * @return          The reference to the current matrix that contains the result.
 */
template< typename eT >
inline
const matrix< complex< eT > >& matrix< complex< eT > >::operator-=(const matrix< complex< eT > >& A)
{
    if ( rows != A.rows || cols != A.cols )
    {
        uzlmath_error("%s", "Dimension mismatch in matrix-matrix subtraction.");
    }
    
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        mem[i] -= A.mem[i];
    }
    
    return *this;
}

/*!
 * @brief           The multiplication assignment operator.
 * @details         Multiplying the matrix \f$A\f$ to the current matrix and stores the
 *                  result in the current matrix. Both matrices should have the same size
 *
 * @param[in]       A The matrix that is supposed to be multiplied to the current matrix
 *
 * @return          The reference to the current matrix that contains the result.
 */
template< typename eT >
inline
const matrix< complex< eT > >& matrix< complex< eT > >::operator*=(const matrix< eT >& A)
{
    if ( cols != A.rows )
    {
        uzlmath_error("%s", "Dimension mismatch in matrix-matrix multiplication.");
    }
    
    size_t n_rows = rows;
    size_t n_cols = A.cols();
    
    if (is_int< eT >::value == true || is_long< eT >::value == true || is_short< eT >::value == true)
    {
        float* tmp_mem  = new float[2 * rows * cols];
        float* tmp_A    = new float[2 * A.rows * A.cols];
        
        size_t i, cap = rows * cols, cap_a = A.rows * A.cols;
        for (i = 0; i < cap; ++i)
        {
            tmp_mem[i * 2]      = static_cast< float >(mem[i].re);
            tmp_mem[i * 2 + 1]  = static_cast< float >(mem[i].im);
        }
        
        for (i = 0; i < cap_a; ++i)
        {
            tmp_A[i * 2]        = static_cast< float >(A.mem[i]);
            tmp_A[i * 2 + 1]    = static_cast< float >(0);
        }
        
        float* C        = new float[2 * n_rows * n_cols];
        float alpha[2]  = {1, 0};
        float beta[2]   = {0, 0};
        uzlblas_cgemm(UZLblasNoTrans, UZLblasNoTrans, n_rows, n_cols, cols, alpha, tmp_mem, rows, tmp_A, A.rows, beta, C, n_rows);
        
        delete [] mem;
        mem = new complex< eT >[n_rows * n_cols];
        
        size_t cap_c = n_rows * n_cols;
        for (i = 0; i < cap_c; ++i)
        {
            mem[i] = complex< eT >(C[i * 2], C[i * 2 + 1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_A;
        delete [] C;
    }
    else if (is_float< eT >::value == true)
    {
        float* tmp_mem  = reinterpret_cast< float* >(mem);
        float* tmp_A    = new float[2 * A.rows * A.cols];
        
        size_t i, cap_a = A.rows * A.cols;
        for (i = 0; i < cap_a; ++i)
        {
            tmp_A[i * 2]        = static_cast< float >(A.mem[i]);
            tmp_A[i * 2 + 1]    = static_cast< float >(0);
        }
        
        float* C        = new float[2 * n_rows * n_cols];
        float alpha[2]  = {1, 0};
        float beta[2]   = {0, 0};
        uzlblas_cgemm(UZLblasNoTrans, UZLblasNoTrans, n_rows, n_cols, cols, alpha, tmp_mem, rows, tmp_A, A.rows, beta, C, n_rows);
        
        delete [] mem;
        mem = new complex< eT >[n_rows * n_cols];
        
        size_t cap_c = 2 * n_rows * n_cols;
        memcpy(mem, C, cap_c * sizeof(float));
        
        delete [] tmp_A;
        delete [] C;
    }
    else if (is_double< eT >::value == true)
    {
        double* tmp_mem = reinterpret_cast< double* >(mem);
        double* tmp_A   = new double[2 * A.rows * A.cols];
        
        size_t i, cap_a = A.rows * A.cols;
        for (i = 0; i < cap_a; ++i)
        {
            tmp_A[i * 2]        = static_cast< double >(A.mem[i]);
            tmp_A[i * 2 + 1]    = static_cast< double >(0);
        }
        
        double* C       = new double[2 * n_rows * n_cols];
        double alpha[2] = {1, 0};
        double beta[2]  = {0, 0};
        uzlblas_zgemm(UZLblasNoTrans, UZLblasNoTrans, n_rows, n_cols, cols, alpha, tmp_mem, rows, tmp_A, A.rows, beta, C, n_rows);
        
        delete [] mem;
        mem = new complex< eT >[n_rows * n_cols];
        
        size_t cap_c = 2 * n_rows * n_cols;
        memcpy(mem, C, cap_c * sizeof(double));
        
        delete [] tmp_A;
        delete [] C;
    }
    else
    {
        double* tmp_mem = new double[2 * rows * cols];
        double* tmp_A   = new double[2 * A.rows * A.cols];
        
        size_t i, cap = rows * cols, cap_a = A.rows * A.cols;
        for (i = 0; i < cap; ++i)
        {
            tmp_mem[i * 2]      = static_cast< double >(mem[i].re);
            tmp_mem[i * 2 + 1]  = static_cast< double >(mem[i].im);
        }
        
        for (i = 0; i < cap_a; ++i)
        {
            tmp_A[i * 2]        = static_cast< double >(A.mem[i]);
            tmp_A[i * 2 + 1]    = static_cast< double >(0);
        }
        
        double* C       = new double[2 * n_rows * n_cols];
        double alpha[2] = {1, 0};
        double beta[2]  = {0, 0};
        uzlblas_zgemm(UZLblasNoTrans, UZLblasNoTrans, n_rows, n_cols, cols, alpha, tmp_mem, rows, tmp_A, A.rows, beta, C, n_rows);
        
        delete [] mem;
        mem = new complex< eT >[n_rows * n_cols];
        
        size_t cap_c = n_rows * n_cols;
        for (i = 0; i < cap_c; ++i)
        {
            mem[i] = complex< eT >(C[i * 2], C[i * 2 + 1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_A;
        delete [] C;
    }
    
    return *this;
}

/*!
 * @brief           The multiplication assignment operator.
 * @details         Multiplying the matrix \f$A\f$ to the current matrix and stores the
 *                  result in the current matrix.
 *
 * @param[in]       A The matrix that is supposed to be multiplied to the current matrix
 *
 * @return          The reference to the current matrix that contains the result.
 */
template< typename eT >
inline
const matrix< complex< eT > >& matrix< complex< eT > >::operator*=(const matrix< complex< eT > >& A)
{
    if ( cols != A.rows )
    {
        uzlmath_error("%s", "Dimension mismatch in matrix-matrix multiplication.");
    }
    
    size_t n_rows = rows;
    size_t n_cols = A.cols();
    
    if (is_int< eT >::value == true || is_long< eT >::value == true || is_short< eT >::value == true)
    {
        float* tmp_mem  = new float[2 * rows * cols];
        float* tmp_A    = new float[2 * A.rows * A.cols];
        
        size_t i, cap = rows * cols, cap_a = A.rows * A.cols;
        for (i = 0; i < cap; ++i)
        {
            tmp_mem[i * 2]      = static_cast< float >(mem[i].re);
            tmp_mem[i * 2 + 1]  = static_cast< float >(mem[i].im);
        }
        
        for (i = 0; i < cap_a; ++i)
        {
            tmp_A[i * 2]        = static_cast< float >(A.mem[i].re);
            tmp_A[i * 2 + 1]    = static_cast< float >(A.mem[i].im);
        }
        
        float* C        = new float[2 * n_rows * n_cols];
        float alpha[2]  = {1, 0};
        float beta[2]   = {0, 0};
        uzlblas_cgemm(UZLblasNoTrans, UZLblasNoTrans, n_rows, n_cols, cols, alpha, tmp_mem, rows, tmp_A, A.rows, beta, C, n_rows);
        
        delete [] mem;
        mem = new complex< eT >[n_rows * n_cols];
        
        size_t cap_c = n_rows * n_cols;
        for (i = 0; i < cap_c; ++i)
        {
            mem[i] = complex< eT >(C[i * 2], C[i * 2 + 1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_A;
        delete [] C;
    }
    if (is_float< eT >::value == true)
    {
        float* tmp_mem = reinterpret_cast< float* >(mem);
        float* tmp_A   = reinterpret_cast< float* >(A.mem);
        
        float* C       = new float[2 * n_rows * n_cols];
        float alpha[2] = {1, 0};
        float beta[2]  = {0, 0};
        uzlblas_cgemm(UZLblasNoTrans, UZLblasNoTrans, n_rows, n_cols, cols, alpha, tmp_mem, rows, tmp_A, A.rows, beta, C, n_rows);
        
        delete [] mem;
        mem = new complex< eT >[n_rows * n_cols];
        
        size_t cap_c = 2 * n_rows * n_cols;
        memcpy(mem, C, cap_c * sizeof(float));
        
        delete [] C;
    }
    if (is_double< eT >::value == true)
    {
        double* tmp_mem = reinterpret_cast< double* >(mem);
        double* tmp_A   = reinterpret_cast< double* >(A.mem);
        
        double* C       = new double[2 * n_rows * n_cols];
        double alpha[2] = {1, 0};
        double beta[2]  = {0, 0};
        uzlblas_zgemm(UZLblasNoTrans, UZLblasNoTrans, n_rows, n_cols, cols, alpha, tmp_mem, rows, tmp_A, A.rows, beta, C, n_rows);
        
        delete [] mem;
        mem = new complex< eT >[n_rows * n_cols];
        
        size_t cap_c = 2 * n_rows * n_cols;
        memcpy(mem, C, cap_c * sizeof(float));
        
        delete [] C;
    }
    else
    {
        double* tmp_mem = new double[2 * rows * cols];
        double* tmp_A   = new double[2 * A.rows * A.cols];
        
        size_t i, cap = rows * cols, cap_a = A.rows * A.cols;
        for (i = 0; i < cap; ++i)
        {
            tmp_mem[i * 2]      = static_cast< double >(mem[i].re);
            tmp_mem[i * 2 + 1]  = static_cast< double >(mem[i].im);
        }
        
        for (i = 0; i < cap_a; ++i)
        {
            tmp_A[i * 2]        = static_cast< double >(A.mem[i].re);
            tmp_A[i * 2 + 1]    = static_cast< double >(A.mem[i].im);
        }
        
        double* C       = new double[2 * n_rows * n_cols];
        double alpha[2] = {1, 0};
        double beta[2]  = {0, 0};
        uzlblas_zgemm(UZLblasNoTrans, UZLblasNoTrans, n_rows, n_cols, cols, alpha, tmp_mem, rows, tmp_A, A.rows, beta, C, n_rows);
        
        delete [] mem;
        mem = new complex< eT >[n_rows * n_cols];
        
        size_t cap_c = n_rows * n_cols;
        for (i = 0; i < cap_c; ++i)
        {
            mem[i] = complex< eT >(C[i * 2], C[i * 2 + 1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_A;
        delete [] C;
    }
    
    return *this;
}



/*!
 * @brief           The multiplication assignment operator for matrix-complex vector product.
 * @details         Multiplying the complex vector \f$v\f$ to the current matrix and stores the
 *                  result in the current matrix.
 *
 * @param[in]       v The complex vector that is supposed to be multiplied to the current matrix
 *
 * @return          The reference to the current matrix that contains the result.
 */
template< typename eT >
inline
const matrix< complex< eT > >& matrix< complex< eT > >::operator*=(const vector< complex< eT > >& v)
{
    if ( v.type == vec_type::ROW || cols != v.n_elements() )
    {
        uzlmath_error("%s", "Dimension mismatch in matrix-matrix multiplication.");
    }
    
    // create new memory array
    complex< eT >* new_mem = new complex< eT >[rows];
    
    // set each value to 0
    memset(new_mem, 0, 2 * rows * sizeof(eT));
    
    // do multiplication
    size_t i, j;
    for (j = 0; j < cols; ++j)
    {
        for (i = 0; i < rows; ++i)
        {
            new_mem[i] +=  mem[j * rows + i] * v[j];
        }
    }
    
    // free current memory array
    delete [] mem;
    
    // adjust size and memory
    cols = 1;
    mem = new_mem;
    
    // return reference to the current matrix
    return *this;
}

/*!
 * @brief           The plus sign operator.
 * @details         Copies the matrix applying
 *                  \f{eqnarray*}{
 *                      B_{ij} = +A_{ij}
 *                  \f}
 *
 * @return          A matrix containing a copy of the current matrix.
 */
template< typename eT >
inline
matrix< complex< eT > > matrix< complex< eT > >::operator+()
{
    matrix< complex< eT > > C(rows, cols);
    memcpy(C.mem, mem, rows * cols * sizeof(complex< eT >));
    
    return C;
}

/*!
 * @brief           The minus sign operator.
 * @details         Copies the matrix applying
 *                  \f{eqnarray*}{
 *                      B_{ij} = -A_{ij}
 *                  \f}
 *
 * @return          A matrix containing the negative values of the current matrix.
 */
template< typename eT >
inline
matrix< complex< eT > > matrix< complex< eT > >::operator-()
{
    matrix< complex< eT > > C(rows, cols);
    
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        C.mem[i] = - mem[i];
    }
    
    return C;
}

/*!
 * @brief           The plus operator for adding a scalar to the matrix.
 * @details         The scalar represents a matrix of the same dimension as this
 *                  matrix and contains only the scalar value on each entry. The
 *                  matrix will be added element-wise to the current matrix.
 *
 * @param[in]       rhs The scalar value
 *
 * @return          A new matrix containing the result of addition.
 *                  \f{eqnarray*}{
 *                      \left(\begin{array}{c c c c}
 *                          a_{11} + rhs & a_{12} + rhs & \cdots & a_{1n} + rhs\\
 *                          a_{21} + rhs & a_{22} + rhs & \cdots & a_{2n} + rhs\\
 *                          \vdots & \vdots & \ddots & \vdots\\
 *                          a_{m1} + rhs & a_{m2} + rhs & \cdots & a_{mn} + rhs
 *                      \end{array}\right)
 *                  \f}
 */
template< typename eT >
inline
matrix< complex< eT > > matrix< complex< eT > >::operator+(const eT& rhs)
{
    matrix< complex< eT > > C(rows, cols);
    
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        C.mem[i] = mem[i] + complex< eT >(rhs, 0);
    }
    
    return C;
}

/*!
 * @brief           The plus operator for adding a scalar to the matrix.
 * @details         The scalar represents a matrix of the same dimension as this
 *                  matrix and contains only the scalar value on each entry. The
 *                  matrix will be added element-wise to the current matrix.
 *
 * @param[in]       rhs The scalar value
 *
 * @return          A new matrix containing the result of addition.
 *                  \f{eqnarray*}{
 *                      \left(\begin{array}{c c c c}
 *                          a_{11} + rhs & a_{12} + rhs & \cdots & a_{1n} + rhs\\
 *                          a_{21} + rhs & a_{22} + rhs & \cdots & a_{2n} + rhs\\
 *                          \vdots & \vdots & \ddots & \vdots\\
 *                          a_{m1} + rhs & a_{m2} + rhs & \cdots & a_{mn} + rhs
 *                      \end{array}\right)
 *                  \f}
 */
template< typename eT >
inline
matrix< complex< eT > > matrix< complex< eT > >::operator+(const complex< eT >& rhs)
{
    matrix< complex< eT > > C(rows, cols);
    
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        C.mem[i] = mem[i] + rhs;
    }
    
    return C;
}

/*!
 * @brief           The minus operator for subtracting a scalar to the matrix.
 * @details         The scalar represents a matrix of the same dimension as this
 *                  matrix and contains only the scalar value on each entry. The
 *                  matrix will be subtracted element-wise to the current matrix.
 *
 * @param[in]       rhs The scalar value
 *
 * @return          A new matrix containing the result of subtraciton.
 *                  \f{eqnarray*}{
 *                      \left(\begin{array}{c c c c}
 *                          a_{11} - rhs & a_{12} - rhs & \cdots & a_{1n} - rhs\\
 *                          a_{21} - rhs & a_{22} - rhs & \cdots & a_{2n} - rhs\\
 *                          \vdots & \vdots & \ddots & \vdots\\
 *                          a_{m1} - rhs & a_{m2} - rhs & \cdots & a_{mn} - rhs
 *                      \end{array}\right)
 *                  \f}
 */
template< typename eT >
inline
matrix< complex< eT > > matrix< complex< eT > >::operator-(const eT& rhs)
{
    matrix< complex< eT > > C(rows, cols);
    
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        C.mem[i] = mem[i] - complex< eT >(rhs, 0);
    }
    
    return C;
}

/*!
 * @brief           The minus operator for subtracting a scalar to the matrix.
 * @details         The scalar represents a matrix of the same dimension as this
 *                  matrix and contains only the scalar value on each entry. The
 *                  matrix will be subtracted element-wise to the current matrix.
 *
 * @param[in]       rhs The scalar value
 *
 * @return          A new matrix containing the result of subtraciton.
 *                  \f{eqnarray*}{
 *                      \left(\begin{array}{c c c c}
 *                          a_{11} - rhs & a_{12} - rhs & \cdots & a_{1n} - rhs\\
 *                          a_{21} - rhs & a_{22} - rhs & \cdots & a_{2n} - rhs\\
 *                          \vdots & \vdots & \ddots & \vdots\\
 *                          a_{m1} - rhs & a_{m2} - rhs & \cdots & a_{mn} - rhs
 *                      \end{array}\right)
 *                  \f}
 */
template< typename eT >
inline
matrix< complex< eT > > matrix< complex< eT > >::operator-(const complex< eT >& rhs)
{
    matrix< complex< eT > > C(rows, cols);
    
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        C.mem[i] = mem[i] - rhs;
    }
    
    return C;
}

/*!
 * @brief           The multiply operator for scaling the matrix elements.
 * @details         The scalar will be multiplied with each element in the current
 *                  matrix.
 *
 * @param[in]       rhs The scalar value
 *
 * @return          A new matrix containing the result of multiplication
 *                  \f{eqnarray*}{
 *                      \lambda \left(\begin{array}{c c c c}
 *                          a_{11} & a_{12} & \cdots & a_{1n}\\
 *                          a_{21} & a_{22} & \cdots & a_{2n}\\
 *                          \vdots & \vdots & \ddots & \vdots\\
 *                          a_{m1} &  a_{m2} & \cdots & a_{mn}
 *                      \end{array}\right) =
 *                      \left(\begin{array}{c c c c}
 *                          \lambda a_{11} & \lambda a_{12} & \cdots & \lambda a_{1n}\\
 *                          \lambda a_{21} & \lambda a_{22} & \cdots & \lambda a_{2n}\\
 *                          \vdots & \vdots & \ddots & \vdots\\
 *                          \lambda a_{m1} & \lambda a_{m2} & \cdots & \lambda a_{mn}
 *                      \end{array}\right)
 *                  \f}
 */
template< typename eT >
inline matrix< complex< eT > > matrix< complex< eT > >::operator*(const eT& rhs)
{
    matrix< complex< eT > > C(rows, cols);
    
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        C.mem[i] = mem[i] * complex< eT >(rhs, 0);
    }
    
    return C;
}

/*!
 * @brief           The multiply operator for scaling the matrix elements.
 * @details         The scalar will be multiplied with each element in the current
 *                  matrix.
 *
 * @param[in]       rhs The scalar value
 *
 * @return          A new matrix containing the result of multiplication
 *                  \f{eqnarray*}{
 *                      \lambda \left(\begin{array}{c c c c}
 *                          a_{11} & a_{12} & \cdots & a_{1n}\\
 *                          a_{21} & a_{22} & \cdots & a_{2n}\\
 *                          \vdots & \vdots & \ddots & \vdots\\
 *                          a_{m1} &  a_{m2} & \cdots & a_{mn}
 *                      \end{array}\right) =
 *                      \left(\begin{array}{c c c c}
 *                          \lambda a_{11} & \lambda a_{12} & \cdots & \lambda a_{1n}\\
 *                          \lambda a_{21} & \lambda a_{22} & \cdots & \lambda a_{2n}\\
 *                          \vdots & \vdots & \ddots & \vdots\\
 *                          \lambda a_{m1} & \lambda a_{m2} & \cdots & \lambda a_{mn}
 *                      \end{array}\right)
 *                  \f}
 */
template< typename eT >
inline matrix< complex< eT > > matrix< complex< eT > >::operator*(const complex< eT >& rhs)
{
    matrix< complex< eT > > C(rows, cols);
    
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        C.mem[i] = mem[i] * rhs;
    }
    
    return C;
}

/*!
 * @brief           The division operator for scaling the matrix elements.
 * @details         The scalar will be divided with each element in the current
 *                  matrix.
 *
 * @param[in]       rhs The scalar value
 *
 * @return          A new matrix containing the result of multiplication
 *                  \f{eqnarray*}{
 *                      \frac{1}{\lambda} \left(\begin{array}{c c c c}
 *                          a_{11} & a_{12} & \cdots & a_{1n}\\
 *                          a_{21} & a_{22} & \cdots & a_{2n}\\
 *                          \vdots & \vdots & \ddots & \vdots\\
 *                          a_{m1} &  a_{m2} & \cdots & a_{mn}
 *                      \end{array}\right) =
 *                      \left(\begin{array}{c c c c}
 *                          \frac{1}{\lambda} a_{11} & \frac{1}{\lambda} a_{12} & \cdots & \frac{1}{\lambda} a_{1n}\\
 *                          \frac{1}{\lambda} a_{21} & \frac{1}{\lambda} a_{22} & \cdots & \frac{1}{\lambda} a_{2n}\\
 *                          \vdots & \vdots & \ddots & \vdots\\
 *                          \frac{1}{\lambda} a_{m1} & \frac{1}{\lambda} a_{m2} & \cdots & \frac{1}{\lambda} a_{mn}
 *                      \end{array}\right)
 *                  \f}
 */
template< typename eT >
inline matrix< complex< eT > > matrix< complex< eT > >::operator/(const eT& rhs)
{
    if ( rhs == 0 )
    {
        uzlmath_error("%s", "Division by zero in matrix-scalar division.");
    }
    
    matrix< complex< eT > > C(rows, cols);
    
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        C.mem[i] = mem[i] / complex< eT >(rhs, 0);
    }
    
    return C;
}

/*!
 * @brief           The division operator for scaling the matrix elements.
 * @details         The scalar will be divided with each element in the current
 *                  matrix.
 *
 * @param[in]       rhs The scalar value
 *
 * @return          A new matrix containing the result of multiplication
 *                  \f{eqnarray*}{
 *                      \frac{1}{\lambda} \left(\begin{array}{c c c c}
 *                          a_{11} & a_{12} & \cdots & a_{1n}\\
 *                          a_{21} & a_{22} & \cdots & a_{2n}\\
 *                          \vdots & \vdots & \ddots & \vdots\\
 *                          a_{m1} &  a_{m2} & \cdots & a_{mn}
 *                      \end{array}\right) =
 *                      \left(\begin{array}{c c c c}
 *                          \frac{1}{\lambda} a_{11} & \frac{1}{\lambda} a_{12} & \cdots & \frac{1}{\lambda} a_{1n}\\
 *                          \frac{1}{\lambda} a_{21} & \frac{1}{\lambda} a_{22} & \cdots & \frac{1}{\lambda} a_{2n}\\
 *                          \vdots & \vdots & \ddots & \vdots\\
 *                          \frac{1}{\lambda} a_{m1} & \frac{1}{\lambda} a_{m2} & \cdots & \frac{1}{\lambda} a_{mn}
 *                      \end{array}\right)
 *                  \f}
 */
template< typename eT >
inline
matrix< complex< eT > > matrix< complex< eT > >::operator/(const complex< eT >& rhs)
{
    if ( rhs == 0 )
    {
        uzlmath_error("%s", "Division by zero in matrix-scalar division.");
    }
    
    matrix< complex< eT > > C(rows, cols);
    
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        C.mem[i] = mem[i] / rhs;
    }
    
    return C;
}

/*!
 * @brief           The power operator for matrices.
 * @details         Uses the multiply operator to calculate the power of a matrix.
 *
 * @param[in]       exp The exponent of the power operation
 *
 * @return          A matrix containing the result of the power operation
 *                  \f{eqnarray*}{
 *                      \left(\begin{array}{c c c c}
 *                          a_{11} & a_{12} & \cdots & a_{1n}\\
 *                          a_{21} & a_{22} & \cdots & a_{2n}\\
 *                          \vdots & \vdots & \ddots & \vdots\\
 *                          a_{m1} &  a_{m2} & \cdots & a_{mn}
 *                      \end{array}\right)^n =
 *                      \underbrace{
 *                      \left(\begin{array}{c c c c}
 *                          a_{11} & a_{12} & \cdots & a_{1n}\\
 *                          a_{21} & a_{22} & \cdots & a_{2n}\\
 *                          \vdots & \vdots & \ddots & \vdots\\
 *                          a_{m1} &  a_{m2} & \cdots & a_{mn}
 *                      \end{array}\right)\cdot\;
 *                      \cdots\;\cdot
 *                      \left(\begin{array}{c c c c}
 *                          a_{11} & a_{12} & \cdots & a_{1n}\\
 *                          a_{21} & a_{22} & \cdots & a_{2n}\\
 *                          \vdots & \vdots & \ddots & \vdots\\
 *                          a_{m1} &  a_{m2} & \cdots & a_{mn}
 *                      \end{array}\right)
 *                      }_{n}
 *                  \f}
 */
template< typename eT >
inline
matrix< complex< eT > > matrix< complex< eT > >::operator^(const unsigned int& exp)
{
    if ( rows != cols )
    {
        uzlmath_error("%s", "Dimension mismatch in matrix power operator.");
    }
    
    matrix< complex< eT > > C = *this;
    unsigned int i = 0;
    for (i = 0; i < exp - 1; ++i)
    {
        C *= *this;
    }
    
    return C;
}



/*!
 * @brief           The addition assignment operator for a scalar value.
 * @details         Adds the scalar value on each entry of the current matrix.
 *
 * @param[in]       rhs The scalar value that is supposed to be added on each entry.
 *
 * @return          The reference to the current matrix.
 */
template< typename eT >
inline
matrix< complex< eT > >& matrix< complex< eT > >::operator+=(const eT& rhs)
{
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        mem[i] += complex< eT >(rhs, 0);
    }
    return *this;
}

/*!
 * @brief           The addition assignment operator for a scalar value.
 * @details         Adds the scalar value on each entry of the current matrix.
 *
 * @param[in]       rhs The scalar value that is supposed to be added on each entry.
 *
 * @return          The reference to the current matrix.
 */
template< typename eT >
inline matrix< complex< eT > >& matrix< complex< eT > >::operator+=(const complex< eT >& rhs)
{
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        mem[i] += rhs;
    }
    return *this;
}

/*!
 * @brief           The addition assignment operator for a scalar value.
 * @details         Adds the scalar value on each entry of the current matrix.
 *
 * @param[in]       rhs The scalar value that is supposed to be added on each entry.
 *
 * @return          The reference to the current matrix.
 */
template< typename eT >
inline
matrix< complex< eT > >& matrix< complex< eT > >::operator-=(const eT& rhs)
{
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        mem[i] -= complex< eT >(rhs, 0);
    }
    return *this;
}

/*!
 * @brief           The addition assignment operator for a scalar value.
 * @details         Adds the scalar value on each entry of the current matrix.
 *
 * @param[in]       rhs The scalar value that is supposed to be added on each entry.
 *
 * @return          The reference to the current matrix.
 */
template< typename eT >
inline matrix< complex< eT > >& matrix< complex< eT > >::operator-=(const complex< eT >& rhs)
{
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        mem[i] -= rhs;
    }
    return *this;
}

/*!
 * @brief           The multiply assignment operator for a scalar value.
 * @details         Multiplies the scalar value to each entry of the current matrix.
 *
 * @param[in]       rhs The scalar value that is supposed to be multiplied to each entry.
 *
 * @return          The reference to the current matrix.
 */
template< typename eT >
inline
matrix< complex< eT > >& matrix< complex< eT > >::operator*=(const eT& rhs)
{
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        mem[i] *= complex< eT >(rhs, 0);
    }
    return *this;
}

/*!
 * @brief           The multiply assignment operator for a scalar value.
 * @details         Multiplies the scalar value to each entry of the current matrix.
 *
 * @param[in]       rhs The scalar value that is supposed to be multiplied to each entry.
 *
 * @return          The reference to the current matrix.
 */
template< typename eT >
inline
matrix< complex< eT > >& matrix< complex< eT > >::operator*=(const complex< eT >& rhs)
{
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        mem[i] *= rhs;
    }
    return *this;
}

/*!
 * @brief           The division assignment operator for a scalar value.
 * @details         Divides each entry of the current matrix by the given scalar value.
 *
 * @param[in]       rhs The scalar value that is used for the division.
 *
 * @return          The reference to the current matrix.
 */
template< typename eT >
inline
matrix< complex< eT > >& matrix< complex< eT > >::operator/=(const eT& rhs)
{
    if ( rhs == 0 )
    {
        uzlmath_error("%s", "Division by zero in matrix scalar division.");
    }
    
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        mem[i] /= complex< eT >(rhs, 0);
    }
    return *this;
}

/*!
 * @brief           The division assignment operator for a scalar value.
 * @details         Divides each entry of the current matrix by the given scalar value.
 *
 * @param[in]       rhs The scalar value that is used for the division.
 *
 * @return          The reference to the current matrix.
 */
template< typename eT >
inline
matrix< complex< eT > >& matrix< complex< eT > >::operator/=(const complex< eT >& rhs)
{
    if ( rhs == 0 )
    {
        uzlmath_error("%s", "Division by zero in matrix scalar division.");
    }
    
    size_t i, cap = rows * cols;
    for (i = 0; i < cap; ++i)
    {
        mem[i] /= rhs;
    }
    return *this;
}

/*!
 * @brief           The power assignment operator for the scalar value.
 * @details         Multiplies the current matrix \f$n\f$ times with itself where \f$n\f$
 *                  is the given parameter _rhs_
 *
 * @param[in]       exp The exponent of the power operation.
 *
 * @return          The reference to the current matrix.
 */
template< typename eT >
inline
matrix< complex< eT > >& matrix< complex< eT > >::operator^=(const unsigned int& exp)
{
    if ( rows != cols )
    {
        uzlmath_error("%s", "Dimension mismatch in matrix power operator.");
    }
    
    matrix< complex< eT > > C = *this;
    unsigned int i;
    for (i = 0; i < exp; ++i)
    {
        *this *= C;
    }
    
    return *this;
}


//inline matrix< complex< eT > >& operator<<(const eT& val);
//inline matrix< complex< eT > >& operator<<(const complex< eT >& val);

/*!
 * @brief           Indexing operator to access the specif matrix element.
 * @details         The indexing operator can be used to set a specific element of the
 *                  current matrix by using two indices \f$m\f$ and \f$n\f$.
 *
 * @param[in]       i The row of the needed element
 * @param[in]       j The column of the needed element
 *
 * @return          The pointer to matrix element in row \f$i\f$ and column \f$j\f$ (\f$a_{ij}\f$).
 */
template< typename eT >
inline
complex< eT >& matrix< complex< eT > >::operator()(const size_t& i, const size_t& j)
{
    return mem[j * rows + i];
}

/*!
 * @brief           Indexing operator to get the specific matrix element.
 * @details         The indexing operator can be used to get a specific element of the
 *                  current matrix by using two indices \f$m\f$ and \f$n\f$.
 *
 * @param[in]       i The row of the needed element
 * @param[in]       j The column of the needed element
 *
 * @return          The matrix element in row \f$i\f$ and column \f$j\f$ (\f$a_{ij}\f$).
 */
template< typename eT >
inline
const complex< eT >& matrix< complex< eT > >::operator()(const size_t& i, const size_t& j) const
{
    return mem[j * rows + i];
}



/*!
 * @brief           Filling the current matrix with zeros.
 * @details         This function will fill the matrix with zeros according to
 *                  \f{eqnarray*}{
 *                      A_{ij} = 0\qquad \forall i,\in \{1,\dots,M\}, j\in\{1,\dots,N\}
 *                  \f}
 *                  where \f$A\f$ is a \f$M\times N\f$ matrix.
 */
template< typename eT >
inline
void matrix< complex< eT > >::zeros()
{
    for (int i = 0; i < rows * cols; ++i)
    {
        mem[i] = complex< eT >(0, 0);
    }
}

/*!
 * @brief           Filling the current matrix with ones.
 * @details         This function will fill the matrix with zeros according to
 *                  \f{eqnarray*}{
 *                      A_{ij} = 1\qquad \forall i,\in \{1,\dots,M\}, j\in\{1,\dots,N\}
 *                  \f}
 *                  where \f$A\f$ is a \f$M\times N\f$ matrix.
 */
template< typename eT >
inline
void matrix< complex< eT > >::ones()
{
    for (int i = 0; i < rows * cols; ++i)
    {
        mem[i] = complex< eT >(1, 0);
    }
}

/*!
 * @brief           Makes the current matrix to a diagonal matrix with ones on its diagonal.
 * @details         This function will fill the matrix with zeros and ones according to
 *                  \f[
 *                      A_{ij} = \begin{cases}
 *                          1 & \text{if } i=j\\
 *                          0 & \text{else }
 *                      \end{cases}
 *                  \f]
 *                  which looks like this
 *                  \f[
 *                      A = \left(\begin{array}{c c c c}
 *                          1 & 0 & \cdots & 0\\
 *                          0 & 1 & \cdots & 0\\
 *                          \vdots & \vdots & \ddots & \vdots\\
 *                          0 & 0 & \cdots & 1
 *                      \end{array}\right)
 *                  \f]
 *                  where \f$A\f$ is a \f$M\times N\f$ matrix.
 */
template< typename eT >
inline
void matrix< complex< eT > >::eye()
{
    zeros();
    
    size_t i, el = (rows < cols) ? rows : cols;
    for (i = 0; i < el; ++i)
    {
        mem[i * rows + i] = complex< eT >(1, 0);
    }
}

/*!
 * @brief           Transposing the current matrix.
 * @details         Swapping the elements of the current matrix according to
 *                  \f[
 *                      A_{ij} = A_{ji}
 *                  \f]
 * @todo            Implement in-place transposing
 */
template< typename eT >
inline
void matrix< complex< eT > >::transpose()
{
    complex< eT >* tmp_mem = new complex< eT >[rows * cols];
    size_t tmp_r = cols;
    size_t tmp_c = rows;
    
    size_t i, j;
    for (i = 0; i < rows; ++i)
    {
        for (j = 0; j < cols; ++j)
        {
            tmp_mem[i * cols + j] = mem[j * rows + i];
        }
    }
    
    delete [] mem;
    rows = tmp_r;
    cols = tmp_c;
    mem = tmp_mem;
}

/*!
 * @brief           Filling the matrix with given scalar value.
 * @details         Each entry of the current matrix gets overwritten with the
 *                  given scalar value.
 */
template< typename eT >
inline
void matrix< complex< eT > >::fill(const eT& value)
{
    size_t i;
    for (i = 0; i < rows * cols; ++i)
    {
        mem[i] = complex< eT >(value, 0);
    }
}

/*!
 * @brief           Filling the matrix with given scalar value.
 * @details         Each entry of the current matrix gets overwritten with the
 *                  given scalar value.
 */
template< typename eT >
inline
void matrix< complex< eT > >::fill(const complex< eT >& value)
{
    std::fill(mem, mem + rows * cols, value);
}

/*!
 * @brief           Getter for the number of rows
 * @details         Returning the number of rows of the current matrix
 *
 * @return          The number of rows
 */
template< typename eT >
inline
constexpr size_t matrix< complex< eT > >::n_rows() const
{
    return rows;
}

/*!
 * @brief           Getter for the number of columns
 * @details         Returning the number of rows of the current matrix
 *
 * @return          The number of columns
 */
template< typename eT >
inline
constexpr size_t matrix< complex< eT > >::n_cols() const
{
    return cols;
}

/*!
 * @brief           Getter for the allocated memory of the current matrix
 * @details         Returning the pointer to the dynamic allocated memory of the
 *                  current matrix
 *
 * @return          The pointer to the allocated memory
 */
template< typename eT >
inline
complex< eT >* matrix< complex< eT > >::memptr()
{
    return mem;
}

/*!
 * @brief           Getter for the allocated memory of the current matrix
 * @details         Returning the pointer to the dynamic allocated memory of the
 *                  current matrix
 *
 * @return          The pointer to the allocated memory
 */
template< typename eT >
inline
const complex< eT >* matrix< complex< eT > >::memptr() const
{
    return mem;
}

/*!
 * @brief           Calculates the determinant of the current matrix
 * @details         Calculating the determinant of a matrix by getting its
 *                  LU-Decomposition. Therefore each matrix can be
 *                  transformed to two matrices that are in lower and upper
 *                  triangular notation such that a given matrix \f$A\f$
 *                  can be written as
 *                  \f[
 *                      P\cdot A = L\cdot U
 *                  \f]
 *                  where \f$P\f$ is a permutation matrix and \f$L\f$ a
 *                  lower triangular matrix and \f$U\f$ an upper triangular
 *                  matrix. The lower triangular matrix only have ones on its
 *                  diagonal elements, so the determinant can be easily
 *                  calculated by multiplying all diagonal elements of the
 *                  \f$U\f$ matrix.
 *
 * @return          The value of the determinant
 */
template< typename eT >
inline
const complex<double> matrix< complex< eT > >::determinant()
{
    if (rows != cols)
    {
        uzlmath_error("%s", "For calculating determinant the matrix must be square.");
    }
    
    size_t i, cap = rows * cols;
    
    double A[2 * cap];
    for (i = 0; i < cap; ++i)
    {
        A[i * 2]        = static_cast< double >(mem[i].re);
        A[i * 2 + 1]    = static_cast< double >(mem[i].im);
    }
    
    /*
     * INFO     = 0:  successful exit
     *          < 0:  if INFO = -i, the i-th argument had an illegal value
     *          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
     *                has been completed, but the factor U is exactly
     *                singular, and division by zero will occur if it is used
     *                to solve a system of equations.
     */
    int INFO = (rows > cols) ? rows : cols;
    int LDA  = rows;
    int IPIV[rows];
    
    int M    = rows;
    int N    = cols;
    
    uzllapack_zgetrf(&M, &N, A, &LDA, IPIV, &INFO);
    
    if (INFO < 0)
    {
        uzlmath_error("LU-Decomposition argument error. Argument number %i was illegal.", -INFO);
    }
    else if (INFO > 0)
    {
        return 0;
    }
    
    double det  = complex< eT >(1, 0);
    int n       = (cols < rows) ? cols : rows;
    for (i = 0; i < n; ++i)
    {
        det *= A[i * rows + i];
    }
    
    return det;
}

/*!
 * @brief           Outstream operator overload for complex matrices.
 * @details         The out-steam operator is used to print the matrix
 *                  in a nice form over the std::cout stream.
 * @param[in,out]   o The stream object
 * @param[in]       A The matrix that should be printed
 *
 * @return          The reference to the given out-stream.
 *
 * @ingroup         matrix
 */
template< typename eT >
std::ostream& operator<<(std::ostream& o, const matrix< complex< eT > >& A)
{
    std::ios::fmtflags f( std::cout.flags() );
    o << std::endl;
    
    int width   = 20;
    auto format = std::fixed;
    
    if (is_float< eT >::value == false && is_double< eT >::value == false && is_ldouble< eT >::value == false)
    {
        width = 10;
    }
    
    // check values
    size_t i, j;
    for (i = 0; i < A.n_rows(); ++i)
    {
        for (j = 0; j < A.n_cols(); ++j)
        {
            complex< eT > c = A(i, j);
            if (std::abs(c.re) >= 10 || std::abs(c.im) >= 10)
            {
                width   = 22;
                format  = std::fixed;
                
                if (is_float< eT >::value == false && is_double< eT >::value == false && is_ldouble< eT >::value == false)
                {
                    width = 12;
                }
            }
            
            if (std::abs(c.re) >= 100 || std::abs(c.im) >= 100)
            {
                width   = 24;
                format  = std::fixed;
                
                if (is_float< eT >::value == false && is_double< eT >::value == false && is_ldouble< eT >::value == false)
                {
                    width = 14;
                }
            }
            
            if (std::abs(c.re) >= 1000 || std::abs(c.im) >= 1000)
            {
                width   = 28;
                format  = std::scientific;
                
                if (is_float< eT >::value == false && is_double< eT >::value == false && is_ldouble< eT >::value == false)
                {
                    width = 18;
                }
            }
        }
    }
    
    // prepare output and print
    for (i = 0; i < A.n_rows(); ++i)
    {
        for (j = 0; j < A.n_cols(); ++j)
        {
            // get entry
            complex< eT > c = A(i, j);
            
            // create string
            std::ostringstream val;
            
            // add real value to string
            val << format << std::setprecision(4) << c.re;
            val << (c.im < 0 ? " - " : " + ") << (c.im == 0 ?  0 : std::abs(c.im)) << "i";
            
            // get string from stream
            std::string str = val.str();
            
            // set filling character
            o << std::setfill(' ') << std::right << std::setw(width) << str;
        }
        
        // line break for next line
        o << std::endl;
    }
    
    std::cout.flags( f );
    return o;
}

UZLMATH_END

#endif
