//
//  cx_vector_def.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 02.06.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_cx_vector_def_hpp
#define UZLMathLib_cx_vector_def_hpp

UZLMATH_BEGIN

template< typename eT >
inline
vector< complex< eT > >::vector()
    : size(0)
    , type(vec_type::ROW)
    , inj(0)
    , mem(nullptr)
{}

template< typename eT >
inline
vector< complex< eT > >::vector(const size_t& s, const vec_type& type)
    : inj(0)
    , size(s)
    , type(type)
{
    mem = new complex< eT >[s];
}

template< typename eT >
inline
vector< complex< eT > >::vector(const size_t& s, const eT& initial, const vec_type& type)
    : inj(0)
    , size(s)
    , type(type)
{
    mem = new complex< eT >[s];
    
    if (size > 0)
    {
        complex< eT >* fill_mem = const_cast< complex< eT >* >(mem);
        std::fill(fill_mem, fill_mem + size, complex< eT >(initial, 0));
    }
}

template< typename eT >
inline
vector< complex< eT > >::vector(const size_t& s, const complex< eT >& initial, const vec_type& type)
    : size(s)
    , type(type)
    , inj(0)
{
    mem = new complex< eT >[s];
    
    if (size > 0)
    {
        complex< eT >* fill_mem = const_cast< complex< eT >* >(mem);
        std::fill(fill_mem, fill_mem + size, initial);
    }
}

template< typename eT >
inline
vector< complex< eT > >::vector(const vector< eT >& vec)
    : size(vec.size)
    , type(vec.type)
    , inj(vec.inj)
{
    mem = new complex< eT >[vec.size];
    
    size_t i;
    for (i = 0; i < vec.size; ++i)
    {
        mem[i] = complex< eT >(vec[i], 0);
    }
}

template< typename eT >
inline
vector< complex< eT > >::vector(const vector< complex< eT > >& vec)
    : size(vec.size)
    , type(vec.type)
    , inj(vec.inj)
{
    mem = new complex< eT >[vec.size];
    memcpy(access::rwp(mem), access::rwp(vec.mem), vec.size * sizeof(complex< eT >));
}

template< typename eT >
inline
vector< complex< eT > >::vector(const vector< eT >& vec, const vec_type& type)
    : size(vec.size)
    , type(type)
    , inj(vec.inj)
{
    mem = new complex< eT >[vec.size];
    
    size_t i;
    for (i = 0; i < vec.size; ++i)
    {
        mem[i] = complex< eT >(vec[i], 0);
    }
}

template< typename eT >
inline
vector< complex< eT > >::vector(const vector< complex< eT > >& vec, const vec_type& t)
    : size(vec.size)
    , type(type)
    , inj(vec.inj)
{
    mem = new complex< eT >[size];
    
    if (size > 0)
    {
        memcpy(access::rw(mem), access::rw(vec.mem), size * sizeof(complex< eT >));
    }
}

template< typename eT >
inline
vector< complex< eT > >::vector(vector< complex< eT > >&& vec)
    : inj(vec.inj)
    , size(vec.size)
    , type(vec.type)
{
    const complex< eT >* tmp = mem;
    mem                      = vec.mem;
    vec.mem                  = tmp;
}

template< typename eT >
inline
vector< complex< eT > >::~vector()
{
    delete [] mem;
}

template< typename eT >
inline
vector< complex< eT > >  vector< complex< eT > >::operator+(const vector< eT >& v)
{
    if ( size != v.size || type != v.type)
    {
        uzlmath_error("%s", "size mismatch in complex vector-vector addition.");
    }
    
    vector< complex< eT > > result(size, type);
    
    size_t i;
    for (i = 0; i < v.size; ++i)
    {
        result[i] = mem[i] + complex< eT >(v[i], 0);
    }
    
    return result;
}

template< typename eT >
inline
vector< complex< eT > > vector< complex< eT > >::operator-(const vector< eT >& v)
{
    if ( size != v.size || type != v.type)
    {
        uzlmath_error("%s", "size mismatch in complex vector-vector subtraction.");
    }
    
    vector< complex< eT > > result(size, type);
    
    size_t i;
    for (i = 0; i < v.size; ++i)
    {
        result[i] = mem[i] + complex< eT >(v[i], 0);
    }
    
    return result;
}

template< typename eT >
inline
matrix< complex< eT > > vector< complex< eT > >::operator*(const vector< eT >& v)
{
    if (type == v.type || (type == vec_type::ROW && size != v.size))
    {
        uzlmath_error("%s", "size mismatch in complex vector-vector multiplication.");
    }
    
    int M   = (type   == vec_type::COLUMN ? size   : 1);
    int N   = (v.type == vec_type::ROW    ? v.size : 1);
    int K   = (type   == vec_type::ROW    ? size   : 1);
    int LDA = M;
    int LDB = K;
    int LDC = M;
    
    matrix< complex< eT > > result(M, N);
    
    // do multiplication
    if (is_int< eT >::value == true || is_short< eT >::value == true)
    {
        float* tmp_mem  = new float[2 * size];
        float* tmp_v    = new float[2 * v.size];
        
        size_t i;
        for (i = 0; i < size; ++i)
        {
            tmp_mem[i * 2]      = static_cast< float >(mem[i].re);
            tmp_mem[i * 2 + 1]  = static_cast< float >(mem[i].im);
        }
        
        for (i = 0; i < v.size; ++i)
        {
            tmp_v[i * 2]        = static_cast< float >(v[i]);
            tmp_v[i * 2 + 1]    = static_cast< float >(0);
        }
        
        float* C        = new float[2 * M * N];
        float alpha[2]  = {1.0, 0.0};
        float beta[2]   = {0.0, 0.0};
        uzl_blas_cgemm(UZLblasNoTrans, UZLblasNoTrans, M, N, K, alpha, tmp_mem, LDA, tmp_v, LDB, beta, C, LDC);
        
        size_t cap_c = M * N;
        for (i = 0; i < cap_c; ++i)
        {
            result.mem[i].re = static_cast< eT >(C[i * 2    ]);
            result.mem[i].im = static_cast< eT >(C[i * 2 + 1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_v;
        delete [] C;
    }
    else
    {
        double* tmp_mem = new double[2 * size];
        double* tmp_v   = new double[2 * v.size];
        
        size_t i;
        for (i = 0; i < size; ++i)
        {
            tmp_mem[i * 2]      = static_cast< double >(mem[i].re);
            tmp_mem[i * 2 + 1]  = static_cast< double >(mem[i].im);
        }
        
        for (i = 0; i < v.size; ++i)
        {
            tmp_v[i * 2]        = static_cast< double >(v[i]);
            tmp_v[i * 2 + 1]    = static_cast< double >(0);
        }
        
        double* C       = new double[2 * M * N];
        double alpha[2] = {1.0, 0.0};
        double beta[2]  = {0.0, 0.0};
        uzl_blas_zgemm(UZLblasNoTrans, UZLblasNoTrans, M, N, K, alpha, tmp_mem, LDA, tmp_v, LDB, beta, C, LDC);
        
        size_t cap_c = M * N;
        for (i = 0; i < cap_c; ++i)
        {
            result.mem[i].re = static_cast< eT >(C[i * 2    ]);
            result.mem[i].im = static_cast< eT >(C[i * 2 + 1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_v;
        delete [] C;
    }
    
    return result;
}

template< typename eT >
inline
vector< complex< eT > > vector< complex< eT > >::operator/(const vector< eT >& v)
{
    if (type != v.type || size != size)
    {
        uzlmath_error("%s", "type or size mismatch in element-wise complex vector division.");
    }
    
    vector< complex< eT > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        if (v[i] == 0)
        {
            uzlmath_error("%s", "division by zero in element-wise complex vector division.");
        }
        
        result[i] = mem[i] / complex< eT >(v[i], 0);
    }
    
    return result;
}

template< typename eT >
inline
vector< complex< eT > > vector< complex< eT > >::operator%(const vector< eT >& v)
{
    if (type != v.type || size != size)
    {
        uzlmath_error("%s", "type or size mismatch in element-wise complex vector multiplication.");
    }
    
    vector< complex< eT > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] * complex< eT >(v[i], 0);
    }
    
    return result;
}

template< typename eT >
inline
vector< complex< eT > > vector< complex< eT > >::operator+(const vector< complex< eT > >& v)
{
    if ( size != v.size || type != v.t)
    {
        uzlmath_error("%s", "size mismatch in complex vector-vector addition.");
    }
    
    vector< complex< eT > > result(size, type);
    
    size_t i;
    for (i = 0; i < v.size; ++i)
    {
        result[i] = mem[i] + v[i];
    }
    
    return result;
}

template< typename eT >
inline
vector< complex< eT > > vector< complex< eT > >::operator-(const vector< complex< eT > >& v)
{
    if ( size != v.size || type != v.type)
    {
        uzlmath_error("%s", "size mismatch in complex vector-vector subtraction.");
    }
    
    vector< complex< eT > > result(size, type);
    
    size_t i;
    for (i = 0; i < v.size; ++i)
    {
        result[i] = mem[i] - v[i];
    }
    
    return result;
}

template< typename eT >
inline
matrix< complex< eT > > vector< complex< eT > >::operator*(const vector< complex< eT > >& v)
{
    if (type == v.type || (type == vec_type::ROW && size != v.size))
    {
        uzlmath_error("%s", "size mismatch in vector-vector multiplication.");
    }
    
    int M   = (type   == vec_type::COLUMN ? size   : 1);
    int N   = (v.type == vec_type::ROW    ? v.size : 1);
    int K   = (type   == vec_type::ROW    ? size   : 1);
    int LDA = M;
    int LDB = K;
    int LDC = M;
    
    matrix< complex< eT > > result(M, N);
    
    // do multiplication
    if (is_int< eT >::value == true || is_short< eT >::value == true || is_float< eT >::value == true)
    {
        
        float* tmp_mem  = new float[2 * size];
        float* tmp_v    = new float[2 * v.size];
        
        size_t i;
        for (i = 0; i < size; ++i)
        {
            tmp_mem[i * 2]      = static_cast< float >(mem[i]);
            tmp_mem[i * 2 + 1]  = static_cast< float >(0);
        }
        
        for (i = 0; i < v.size; ++i)
        {
            tmp_v[i * 2]        = static_cast< float >(v[i].re);
            tmp_v[i * 2 + 1]    = static_cast< float >(v[i].im);
        }
        
        float* C        = new float[2 * M * N];
        float alpha[2]  = {1.0, 0.0};
        float beta[2]   = {0.0, 0.0};
        uzl_blas_cgemm(UZLblasNoTrans, UZLblasNoTrans, M, N, K, alpha, tmp_mem, LDA, tmp_v, LDB, beta, C, LDC);
        
        size_t cap_c = M * N;
        for (i = 0; i < cap_c; ++i)
        {
            result.mem[i].re = static_cast< eT >(C[i * 2]);
            result.mem[i].im = static_cast< eT >(C[i * 2 + 1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_v;
        delete [] C;
    }
    else if (is_float< eT >::value == true)
    {
        float alpha[2]   = {1, 0};
        float beta[2]    = {0, 0};
        
        float* A_mem_ptr = reinterpret_cast< float* >(mem);
        float* B_mem_ptr = reinterpret_cast< float* >(v.mem);
        float* C_mem_ptr = reinterpret_cast< float* >(result.mem);
        
        uzl_blas_cgemm(UZLblasNoTrans, UZLblasNoTrans, M, N, K, alpha, A_mem_ptr, LDA, B_mem_ptr, LDB, beta, C_mem_ptr, LDC);
    }
    else if (is_double< eT >::value == true)
    {
        double alpha[2]   = {1, 0};
        double beta[2]    = {0, 0};
        
        double* A_mem_ptr = reinterpret_cast< double* >(mem);
        double* B_mem_ptr = reinterpret_cast< double* >(v.mem);
        double* C_mem_ptr = reinterpret_cast< double* >(result.mem);
        
        uzl_blas_zgemm(UZLblasNoTrans, UZLblasNoTrans, M, N, K, alpha, A_mem_ptr, LDA, B_mem_ptr, LDB, beta, C_mem_ptr, LDC);
    }
    else
    {
        double* tmp_mem = new double[2 * size];
        double* tmp_v   = new double[2 * v.size];
        
        size_t i;
        for (i = 0; i < size; ++i)
        {
            tmp_mem[i * 2]      = static_cast< double >(mem[i]);
            tmp_mem[i * 2 + 1]  = static_cast< double >(0);
        }
        
        for (i = 0; i < v.size; ++i)
        {
            tmp_v[i * 2]        = static_cast< double >(v[i].re);
            tmp_v[i * 2 + 1]    = static_cast< double >(v[i].im);
        }
        
        double* C       = new double[2 * M * N];
        double alpha[2] = {1.0, 0.0};
        double beta[2]  = {0.0, 0.0};
        uzl_blas_zgemm(UZLblasNoTrans, UZLblasNoTrans, M, N, K, alpha, tmp_mem, LDA, tmp_v, LDB, beta, C, LDC);
        
        size_t cap_c = M * N;
        for (i = 0; i < cap_c; ++i)
        {
            result.mem[i].re = static_cast< eT >(C[i * 2]);
            result.mem[i].im = static_cast< eT >(C[i * 2 + 1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_v;
        delete [] C;
    }
    
    return result;
}

template< typename eT >
inline
vector< complex< eT > > vector< complex< eT > >::operator/(const vector< complex< eT > >& v)
{
    if (type != v.type || size != size)
    {
        uzlmath_error("%s", "type or size mismatch in element-wise complex vector division.");
    }
    
    vector< complex< eT > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        if (v[i] == 0)
        {
            uzlmath_error("%s", "division by zero in element-wise complex vector division.");
        }
        
        result[i] = mem[i] / v[i];
    }
    
    return result;
}

template< typename eT >
inline
vector< complex< eT > > vector< complex< eT > >::operator%(const vector< complex< eT > >& v)
{
    if (type != v.type || size != size)
    {
        uzlmath_error("%s", "type or size mismatch in element-wise complex vector multiplication.");
    }
    
    vector< complex< eT > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] * v[i];
    }
    
    return result;
}

template< typename eT >
inline
vector< complex< eT > > vector< complex< eT > >::operator+(const eT& s)
{
    vector< complex< eT > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] + complex< eT >(s, 0);
    }
    
    return result;
}

template< typename eT >
inline
vector< complex< eT > > vector< complex< eT > >::operator-(const eT& s)
{
    vector< complex< eT > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] - complex< eT >(s, 0);
    }
    
    return result;
}

template< typename eT >
inline
vector< complex< eT > > vector< complex< eT > >::operator*(const eT& s)
{
    vector< complex< eT >>  result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = const_cast< complex< double >& >(mem[i]) * complex< eT >(s, 0);
    }
    
    return result;
}

template< typename eT >
inline
vector< complex< eT > > vector< complex< eT > >::operator/(const eT& s)
{
    if (s == 0)
    {
        uzlmath_error("%s", "division by zero in vector-scalar division.");
    }
    
    vector< complex< eT > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] / complex< eT >(s, 0);
    }
    
    return result;
}

template< typename eT >
inline
vector< complex< eT > > vector< complex< eT > >::operator+(const complex< eT >& s)
{
    vector< complex< eT > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] + s;
    }
    
    return result;
}

template< typename eT >
inline
vector< complex< eT > > vector< complex< eT > >::operator-(const complex< eT >& s)
{
    vector< complex< eT > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] - s;
    }
    
    return result;
}

template< typename eT >
inline
vector< complex< eT > > vector< complex< eT > >::operator*(const complex< eT >& s)
{
    vector< complex< eT > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] * s;
    }
    
    return result;
}

template< typename eT >
inline
vector< complex< eT > > vector< complex< eT > >::operator/(const complex< eT >& s)
{
    if (s.re == 0 && s.im == 0)
    {
        uzlmath_error("%s", "division by zero in vector-scalar division.");
    }
    
    vector< complex< eT > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] / s;
    }
    
    return result;
}

template< typename eT >
inline
vector< complex< eT > > vector< complex< eT > >::operator+()
{
    vector< complex< eT > > result(size, type);
    memcpy(access::rw(mem), access::rw(result.mem), size * sizeof(complex< eT >));
    
    return result;
}

template< typename eT >
inline
vector< complex< eT > > vector< complex< eT > >::operator-()
{
    vector< complex< eT > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = - mem[i];
    }
    
    return result;
}


template< typename eT >
inline
vector< complex< eT > > vector< complex< eT > >::operator*(const matrix< eT >& mat)
{
    if ((type == vec_type::ROW && size != mat.rows) || (type == vec_type::COLUMN && mat.rows != 1))
    {
        uzlmath_error("%s", "size mismatch in vector-matrix multiplication.");
    }
    
    vector< complex< eT > > result(mat.n_cols(), vec_type::ROW);
    
    int M   = (type == vec_type::ROW ? 1 : size);
    int N   = mat.n_cols();
    int K   = (type == vec_type::ROW ? size : 1);
    int LDA = M;
    int LDB = K;
    int LDC = M;
    
    if (is_int< eT >::value == true || is_short< eT >::value == true)
    {
        float* tmp_mem  = new float[2 * size];
        float* tmp_mat  = new float[2 * mat.n_cols() * mat.n_rows()];
        
        size_t i;
        for (i = 0; i < size; ++i)
        {
            tmp_mem[i * 2]      = static_cast< float >(mem[i].re);
            tmp_mem[i * 2 + 1]  = static_cast< float >(mem[i].im);
        }
        
        for (i = 0; i < mat.n_cols() * mat.n_rows(); ++i)
        {
            tmp_mat[i * 2]      = static_cast< float >(mat.mem[i]);
            tmp_mat[i * 2 + 1]  = static_cast< float >(0);
        }
        
        float alpha[2]  = {1, 0};
        float beta[2]   = {0, 0};
        float* C        = new float[2 * M * N];
        
        uzl_blas_cgemm(UZLblasNoTrans, UZLblasNoTrans, M, N, K, alpha, tmp_mem, LDA, tmp_mem, LDB, beta, C, LDC);
        
        size_t cap_c = M * N;
        for (i = 0; i < cap_c; ++i)
        {
            result[i].re = static_cast< eT >(C[i * 2    ]);
            result[i].im = static_cast< eT >(C[i * 2 + 1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_mat;
        delete [] C;
    }
    else
    {
        double* tmp_mem = new double[2 * size];
        double* tmp_mat = new double[2 * mat.n_cols() * mat.n_rows()];
        
        size_t i;
        for (i = 0; i < size; ++i)
        {
            tmp_mem[i * 2]      = static_cast< double >(mem[i].re);
            tmp_mem[i * 2 + 1]  = static_cast< double >(mem[i].im);
        }
        
        for (i = 0; i < mat.n_cols() * mat.n_rows(); ++i)
        {
            tmp_mat[i * 2]      = static_cast< double >(mat.mem[i]);
            tmp_mat[i * 2 + 1]  = static_cast< double >(0);
        }
        
        double alpha[2] = {1, 0};
        double beta[2]  = {0, 0};
        double* C       = new double[2 * M * N];
        
        uzl_blas_zgemm(UZLblasNoTrans, UZLblasNoTrans, M, N, K, alpha, tmp_mem, LDA, tmp_mem, LDB, beta, C, LDC);
        
        size_t cap_c = M * N;
        for (i = 0; i < cap_c; ++i)
        {
            result[i].re = static_cast< eT >(C[i * 2    ]);
            result[i].im = static_cast< eT >(C[i * 2 + 1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_mat;
        delete [] C;
    }
    
    return result;
}

template< typename eT >
inline
vector< complex< eT > > vector< complex< eT > >::operator*(const matrix< complex< eT > >& mat)
{
    if ((type == vec_type::ROW && size != mat.rows) || (type == vec_type::COLUMN && mat.rows != 1))
    {
        uzlmath_error("%s", "size mismatch in vector-matrix multiplication.");
    }
    
    vector< complex< eT > > result(mat.n_cols(), vec_type::ROW);
    
    int M   = (type == vec_type::ROW ? 1 : size);
    int N   = mat.n_cols();
    int K   = (type == vec_type::ROW ? size : 1);
    int LDA = M;
    int LDB = K;
    int LDC = M;
    
    if (is_int< eT >::value == true || is_short< eT >::value == true)
    {
        float* tmp_mem  = new float[2 * size];
        float* tmp_mat  = new float[2 * mat.n_cols() * mat.n_rows()];
        
        size_t i;
        for (i = 0; i < size; ++i)
        {
            tmp_mem[i * 2]      = static_cast< float >(mem[i].re);
            tmp_mem[i * 2 + 1]  = static_cast< float >(mem[i].im);
        }
        
        for (i = 0; i < mat.n_cols() * mat.n_rows(); ++i)
        {
            tmp_mat[i * 2]      = static_cast< float >(mat.mem[i].re);
            tmp_mat[i * 2 + 1]  = static_cast< float >(mat.mem[i].im);
        }
        
        float alpha[2]  = {1, 0};
        float beta[2]   = {0, 0};
        float* C        = new float[2 * M * N];
        
        uzl_blas_cgemm(UZLblasNoTrans, UZLblasNoTrans, M, N, K, alpha, tmp_mem, LDA, tmp_mem, LDB, beta, C, LDC);
        
        size_t cap_c = M * N;
        for (i = 0; i < cap_c; ++i)
        {
            result[i].re = static_cast< eT >(C[i * 2    ]);
            result[i].im = static_cast< eT >(C[i * 2 + 1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_mat;
        delete [] C;
    }
    else
    {
        double* tmp_mem = new double[2 * size];
        double* tmp_mat = new double[2 * mat.n_cols() * mat.n_rows()];
        
        size_t i;
        for (i = 0; i < size; ++i)
        {
            tmp_mem[i * 2]      = static_cast< double >(mem[i].re);
            tmp_mem[i * 2 + 1]  = static_cast< double >(mem[i].im);
        }
        
        for (i = 0; i < mat.n_cols() * mat.n_rows(); ++i)
        {
            tmp_mat[i * 2]      = static_cast< double >(mat.mem[i].re);
            tmp_mat[i * 2 + 1]  = static_cast< double >(mat.mem[i].im);
        }
        
        double alpha[2] = {1, 0};
        double beta[2]  = {0, 0};
        double* C       = new double[2 * M * N];
        
        uzl_blas_zgemm(UZLblasNoTrans, UZLblasNoTrans, M, N, K, alpha, tmp_mem, LDA, tmp_mem, LDB, beta, C, LDC);
        
        size_t cap_c = M * N;
        for (i = 0; i < cap_c; ++i)
        {
            result[i].re = static_cast< eT >(C[i * 2    ]);
            result[i].im = static_cast< eT >(C[i * 2 + 1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_mat;
        delete [] C;
    }
    
    return result;
}

template< typename eT >
inline
const vector< complex< eT > >& vector< complex< eT > >::operator=(const vector< eT >& v)
{
    size = v.size;
    type = v.type;
    
    delete [] mem;
    mem = new complex< eT >[size];
    
    if (size > 0)
    {
        size_t i;
        for (i = 0; i < size; ++i)
        {
            mem[i] = complex< eT >(v[i], 0);
        }
    }
    
    return *this;
}

template< typename eT >
inline
const vector< complex< eT > >& vector< complex< eT > >::operator=(const vector< complex< eT > >& v)
{
    if ( this == &v )
    {
        return *this;
    }
    
    size = v.size;
    type = v.type;
    
    delete [] mem;
    mem = new complex< eT >[size];
    
    if (size > 0)
    {
        memcpy(access::rw(mem), access::rw(v.mem), size * sizeof(complex< eT >));
    }
    
    return *this;
}

template< typename eT >
inline
const vector< complex< eT > >& vector< complex< eT > >::operator=(vector< complex< eT > >&& v)
{
    if ( this == &v )
    {
        return *this;
    }
    
    access::rw(size) = v.size;
    access::rw(type) = v.type;
    
    const complex< eT >* tmp = mem;
    mem                      = v.mem;
    v.mem                    = tmp;
    
    return *this;
}


template< typename eT >
inline
const vector< complex< eT > >& vector< complex< eT > >::operator+=(const vector< eT >& v)
{
    if (size != size || type != v.type)
    {
        uzlmath_error("%s", "dimension or size mismatch in complex vector-vector multiplication.");
    }
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] += complex< eT >(v[i], 0);
    }
    
    return *this;
}

template< typename eT >
inline
const vector< complex< eT > >& vector< complex< eT > >::operator+=(const vector< complex< eT > >& v)
{
    if (size != size || type != v.type)
    {
        uzlmath_error("%s", "dimension or size mismatch in complex vector-vector multiplication.");
    }
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] += v[i];
    }
    
    return *this;
}

template< typename eT >
inline
const vector< complex< eT > >& vector< complex< eT > >::operator-=(const vector< eT >& v)
{
    if (size != size || type != v.type)
    {
        uzlmath_error("%s", "dimension or size mismatch in complex vector-vector multiplication.");
    }
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] -= complex< eT >(v[i], 0);
    }
    
    return *this;
}

template< typename eT >
inline
const vector< complex< eT > >& vector< complex< eT > >::operator-=(const vector< complex< eT > >& v)
{
    if (size != size || type != v.type)
    {
        uzlmath_error("%s", "dimension or size mismatch in complex vector-vector multiplication.");
    }
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] -= v[i];
    }
    
    return *this;
}

template< typename eT >
inline
const vector< complex< eT > >& vector< complex< eT > >::operator/=(const vector< eT >& v)
{
    if (type != v.type || size != size)
    {
        uzlmath_error("%s", "type or size mismatch in element-wise complex vector division.");
    }
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        if (v[i] == 0)
        {
            uzlmath_error("%s", "division by zero in element-wise complex vector division.");
        }
        
        mem[i] /= complex< eT >(v[i], 0);
    }
    
    return *this;
}

template< typename eT >
inline
const vector< complex< eT > >& vector< complex< eT > >::operator/=(const vector< complex< eT > >& v)
{
    if (type != v.type || size != size)
    {
        uzlmath_error("%s", "type or size mismatch in element-wise complex vector division.");
    }
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        if (v[i] == 0)
        {
            uzlmath_error("%s", "division by zero in element-wise complex vector division.");
        }
        
        mem[i] /= v[i];
    }
    
    return *this;
}

template< typename eT >
inline
const vector< complex< eT > >& vector< complex< eT > >::operator%=(const vector< eT >& v)
{
    if (type != v.type || size != size)
    {
        uzlmath_error("%s", "type or size mismatch in element-wise complex vector multiplication.");
    }
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] *= complex< eT >(v[i], 0);
    }
    
    return *this;
}

template< typename eT >
inline
const vector< complex< eT > >& vector< complex< eT > >::operator%=(const vector< complex< eT > >& v)
{
    if (type != v.type || size != size)
    {
        uzlmath_error("%s", "type or size mismatch in element-wise complex vector multiplication.");
    }
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] *= v[i];
    }
    
    return *this;
}

template< typename eT >
inline
const vector< complex< eT > >& vector< complex< eT > >::operator*=(const vector< eT >& v)
{
    if (type == v.type || type == vec_type::COLUMN || (type == vec_type::ROW && size != v.size))
    {
        uzlmath_error("%s", "size mismatch or wrong vector type in complex vector-vector multiplication.");
    }
    
    size_t i;
    eT sum = 0;
    for (i = 0; i < size; ++i)
    {
        sum += mem[i] * complex< eT >(v[i], 0);
    }
    
    delete [] mem;
    mem = new complex< eT >[1];
    mem[0] = sum;
    
    return *this;
}

template< typename eT >
inline
const vector< complex< eT > >& vector< complex< eT > >::operator*=(const vector< complex< eT > >& v)
{
    if (type == v.type || type == vec_type::COLUMN || (type == vec_type::ROW && size != v.size))
    {
        uzlmath_error("%s", "size mismatch or wrong vector type in complex vector-vector multiplication.");
    }
    
    size_t i;
    eT sum = 0;
    for (i = 0; i < size; ++i)
    {
        sum += mem[i] * v[i];
    }
    
    delete [] mem;
    mem = new complex< eT >[1];
    mem[0] = sum;
    
    return *this;
}

template< typename eT >
inline
const vector< complex< eT > >& vector< complex< eT > >::operator+=(const eT& s)
{
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] += complex< eT >(s, 0);
    }
    
    return *this;
}

template< typename eT >
inline
const vector< complex< eT > >& vector< complex< eT > >::operator-=(const eT& s)
{
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] -= complex< eT >(s, 0);
    }
    
    return *this;
}

template< typename eT >
inline
const vector< complex< eT > >& vector< complex< eT > >::operator*=(const eT& s)
{
    size_t i;
    for (i = 0; i < size; ++i)
    {
        access::rw(mem[i]) *= complex< eT >(s, 0);
    }
    
    return *this;
}

template< typename eT >
inline
const vector< complex< eT > >& vector< complex< eT > >::operator/=(const eT& s)
{
    if (s == 0)
    {
        uzlmath_error("%s", "division by zero in complex vector-scalar division.");
    }
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] /= complex< eT >(s, 0);
    }
    
    return *this;
}

template< typename eT >
inline
bool vector< complex< eT > >::operator==(const vector< eT >& v)
{
    bool equal = true;
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        if (mem[i] != complex< eT >(v[i], 0))
        {
            equal = false;
            break;
        }
    }
    
    return equal;
}

template< typename eT >
inline
bool vector< complex< eT > >::operator!=(const vector< eT >& v)
{
    return !(*this == v);
}

template< typename eT >
inline
bool vector< complex< eT > >::operator==(const vector< complex< eT > >& v)
{
    bool equal = true;
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        if (mem[i] != v[i])
        {
            equal = false;
            break;
        }
    }
    
    return equal;
}

template< typename eT >
inline
bool vector< complex< eT > >::operator!=(const vector< complex< eT > >& v)
{
    return !(*this == v);
}

template< typename eT >
inline
complex< eT >& vector< complex< eT > >::operator[](const size_t& idx)
{
    return access::rw(mem[idx]);
}

template< typename eT >
inline
constexpr complex< eT >& vector< complex< eT > >::operator[](const size_t& idx) const
{
    return access::rw(mem[idx]);
}

template< typename eT >
inline
complex< eT >& vector< complex< eT > >::operator()(const size_t& idx)
{
    return mem[idx];
}

template< typename eT >
inline
constexpr complex< eT >& vector< complex< eT > >::operator()(const size_t& idx) const
{
    return mem[idx];
}

template< typename eT >
inline
void vector< complex< eT > >::ones()
{
    complex< eT >* fill_mem = const_cast< complex< eT >* >(mem);
    std::fill(fill_mem, fill_mem + size, complex< eT >(1, 0));
}

template< typename eT >
inline
void vector< complex< eT > >::zeros()
{
    complex< eT >* fill_mem = const_cast< complex< eT >* >(mem);
    std::fill(fill_mem, fill_mem + size, complex< eT >(0, 0));
}

template< typename eT >
inline
void vector< complex< eT > >::transpose()
{
    type = (type == vec_type::ROW ? vec_type::COLUMN : vec_type::ROW);
}

template< typename eT >
inline
void vector< complex< eT > >::fill(const eT& s)
{
    complex< eT >* fill_mem = const_cast< complex< eT >* >(mem);
    std::fill(fill_mem, fill_mem + size, complex< eT >(s, 0));
}



template< typename S >
std::ostream& operator<<(std::ostream& o, const vector< complex< S > >& v)
{
    // setting decimal precesion
    std::ios::fmtflags f( std::cout.flags() );
    o << std::endl;
    
    int width   = 20;
    auto format = std::fixed;
    
    if (is_float< S >::value == false && is_double< S >::value == false && is_ldouble< S >::value == false)
    {
        width = 10;
    }
    
    // check values
    size_t i;
    for (i = 0; i < v.n_elements(); ++i)
    {
        complex< S > c = v[i];
        
        if (UZL_ABS(c.re) >= 10 || UZL_ABS(c.im) >= 10)
        {
            width   = 22;
            format  = std::fixed;
            
            if (is_float< S >::value == false && is_double< S >::value == false && is_ldouble< S >::value == false)
            {
                width = 12;
            }
        }
        
        if (UZL_ABS(c.re) >= 100 || UZL_ABS(c.im) >= 100)
        {
            width   = 24;
            format  = std::fixed;
            
            if (is_float< S >::value == false && is_double< S >::value == false && is_ldouble< S >::value == false)
            {
                width = 14;
            }
        }
        
        if (UZL_ABS(c.re) >= 1000 || UZL_ABS(c.im) >= 1000)
        {
            width   = 28;
            format  = std::scientific;
            
            if (is_float< S >::value == false && is_double< S >::value == false && is_ldouble< S >::value == false)
            {
                width = 18;
            }
        }
    }
    
    // prepare output and print
    if (v.vecType() == vec_type::ROW)
    {
        for (i = 0; i < v.n_elements(); ++i)
        {
            // get entry
            complex< S > c = v[i];
            
            // create string
            std::ostringstream val;
            
            // add real value to string
            val << format << std::setprecision(4) << c.re;
            val << (c.im < 0 ? " - " : " + ") << (c.im == 0 ?  0 : UZL_ABS(c.im)) << "i";
            
            // get string from stream
            std::string str = val.str();
            
            // set filling character
            o << std::setfill(' ') << std::right << std::setw(width) << str;
        }
        o << std::endl;
    }
    else
    {
        for (i = 0; i < v.n_elements(); ++i)
        {
            // get entry
            complex< S > c = v[i];
            
            // create string
            std::ostringstream val;
            
            // add real value to string
            val << format << std::setprecision(4) << c.re;
            val << (c.im < 0 ? " - " : " + ") << (c.im == 0 ?  0 : UZL_ABS(c.im)) << "i";
            
            // get string from stream
            std::string str = val.str();
            
            // set filling character
            o << std::setfill(' ') << std::right << std::setw(width) << str << std::endl;
        }
    }
    
    std::cout.flags( f );
    return o;
}

UZLMATH_END

#endif
