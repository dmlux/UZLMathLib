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

template< typename T >
inline
vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::vector()
    : size(0)
    , type(vec_type::ROW)
    , inj(0)
    , mem(nullptr)
{}

template< typename T >
inline
vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::vector(const size_t& s, const vec_type& type)
    : inj(0)
    , size(s)
    , type(type)
{
    mem = new complex< T >[s];
}

template< typename T >
inline
vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::vector(const size_t& s, const T& initial, const vec_type& type)
    : inj(0)
    , size(s)
    , type(type)
{
    mem = new complex< T >[s];
    
    if (size > 0)
    {
        complex< T >* fill_mem = const_cast< complex< T >* >(mem);
        std::fill(fill_mem, fill_mem + size, complex< T >(initial, 0));
    }
}

template< typename T >
inline
vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::vector(const size_t& s, const complex< T >& initial, const vec_type& type)
    : size(s)
    , type(type)
    , inj(0)
{
    mem = new complex< T >[s];
    
    if (size > 0)
    {
        complex< T >* fill_mem = const_cast< complex< T >* >(mem);
        std::fill(fill_mem, fill_mem + size, initial);
    }
}

template< typename T >
inline
vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::vector(const vector< T >& vec)
    : size(vec.size)
    , type(vec.type)
    , inj(vec.inj)
{
    mem = new complex< T >[vec.size];
    
    size_t i;
    for (i = 0; i < vec.size; ++i)
    {
        mem[i] = complex< T >(vec[i], 0);
    }
}

template< typename T >
inline
vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::vector(const vector< complex< T > >& vec)
    : inj(vec.inj)
    , size(vec.size)
    , type(vec.type)
{
    mem = new complex< T >[vec.size];
    memcpy(access::rwp(mem), access::rwp(vec.mem), vec.size * sizeof(complex< T >));
}

template< typename T >
inline
vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::vector(const vector< T >& vec, const vec_type& type)
    : size(vec.size)
    , type(type)
    , inj(vec.inj)
{
    mem = new complex< T >[vec.size];
    
    size_t i;
    for (i = 0; i < vec.size; ++i)
    {
        mem[i] = complex< T >(vec[i], 0);
    }
}

template< typename T >
inline
vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::vector(const vector< complex< T > >& vec, const vec_type& t)
    : size(vec.size)
    , type(type)
    , inj(vec.inj)
{
    mem = new complex< T >[size];
    
    if (size > 0)
    {
        memcpy(access::rw(mem), access::rw(vec.mem), size * sizeof(complex< T >));
    }
}

template< typename T >
inline
vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::vector(vector< complex< T > >&& vec)
    : inj(vec.inj)
    , size(vec.size)
    , type(vec.type)
{
    const complex< T >* tmp = mem;
    mem                      = vec.mem;
    vec.mem                  = tmp;
}

template< typename T >
inline
vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::~vector()
{
    delete [] mem;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator+(const vector< T >& v)
{
    if ( size != v.size || type != v.type)
    {
        uzlmath_error("%s", "size mismatch in complex vector-vector addition.");
    }
    
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < v.size; ++i)
    {
        result[i] = mem[i] + complex< T >(v[i], 0);
    }
    
    return result;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator-(const vector< T >& v)
{
    if ( size != v.size || type != v.type)
    {
        uzlmath_error("%s", "size mismatch in complex vector-vector subtraction.");
    }
    
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < v.size; ++i)
    {
        result[i] = mem[i] + complex< T >(v[i], 0);
    }
    
    return result;
}

template< typename T >
inline
matrix< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator*(const vector< T >& v)
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
    
    matrix< complex< T > > result(M, N);
    
    // do multiplication
    if ( same_type< T, int >::value || same_type< T, short >::value )
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
            result.mem[i].re = static_cast< T >(C[i * 2    ]);
            result.mem[i].im = static_cast< T >(C[i * 2 + 1]);
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
            result.mem[i].re = static_cast< T >(C[i * 2    ]);
            result.mem[i].im = static_cast< T >(C[i * 2 + 1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_v;
        delete [] C;
    }
    
    return result;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator/(const vector< T >& v)
{
    if (type != v.type || size != size)
    {
        uzlmath_error("%s", "type or size mismatch in element-wise complex vector division.");
    }
    
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        if (v[i] == 0)
        {
            uzlmath_error("%s", "division by zero in element-wise complex vector division.");
        }
        
        result[i] = mem[i] / complex< T >(v[i], 0);
    }
    
    return result;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator%(const vector< T >& v)
{
    if (type != v.type || size != size)
    {
        uzlmath_error("%s", "type or size mismatch in element-wise complex vector multiplication.");
    }
    
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] * complex< T >(v[i], 0);
    }
    
    return result;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator+(const vector< complex< T > >& v)
{
    if ( size != v.size || type != v.t)
    {
        uzlmath_error("%s", "size mismatch in complex vector-vector addition.");
    }
    
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < v.size; ++i)
    {
        result[i] = mem[i] + v[i];
    }
    
    return result;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator-(const vector< complex< T > >& v)
{
    if ( size != v.size || type != v.type)
    {
        uzlmath_error("%s", "size mismatch in complex vector-vector subtraction.");
    }
    
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < v.size; ++i)
    {
        result[i] = v[i] - mem[i];
    }
    
    return result;
}

template< typename T >
inline
matrix< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator*(const vector< complex< T > >& v)
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
    
    matrix< complex< T > > result(M, N);
    
    // do multiplication
    if ( same_type< T, int >::value || same_type< T, short >::value || same_type< T, float >::value )
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
            result.mem[i].re = static_cast< T >(C[i * 2]);
            result.mem[i].im = static_cast< T >(C[i * 2 + 1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_v;
        delete [] C;
    }
    else if ( same_type< T, float >::value )
    {
        float alpha[2]   = {1, 0};
        float beta[2]    = {0, 0};
        
        float* A_mem_ptr = reinterpret_cast< float* >(mem);
        float* B_mem_ptr = reinterpret_cast< float* >(v.mem);
        float* C_mem_ptr = reinterpret_cast< float* >(result.mem);
        
        uzl_blas_cgemm(UZLblasNoTrans, UZLblasNoTrans, M, N, K, alpha, A_mem_ptr, LDA, B_mem_ptr, LDB, beta, C_mem_ptr, LDC);
    }
    else if ( same_type< T, double >::value )
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
            result.mem[i].re = static_cast< T >(C[i * 2]);
            result.mem[i].im = static_cast< T >(C[i * 2 + 1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_v;
        delete [] C;
    }
    
    return result;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator/(const vector< complex< T > >& v)
{
    if (type != v.type || size != size)
    {
        uzlmath_error("%s", "type or size mismatch in element-wise complex vector division.");
    }
    
    vector< complex< T > > result(size, type);
    
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

template< typename T >
inline
vector< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator%(const vector< complex< T > >& v)
{
    if (type != v.type || size != size)
    {
        uzlmath_error("%s", "type or size mismatch in element-wise complex vector multiplication.");
    }
    
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] * v[i];
    }
    
    return result;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator+(const T& s)
{
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] + complex< T >(s, 0);
    }
    
    return result;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator-(const T& s)
{
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] - complex< T >(s, 0);
    }
    
    return result;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator*(const T& s)
{
    vector< complex< T >>  result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = const_cast< complex< double >& >(mem[i]) * complex< T >(s, 0);
    }
    
    return result;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator/(const T& s)
{
    if (s == 0)
    {
        uzlmath_error("%s", "division by zero in vector-scalar division.");
    }
    
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] / complex< T >(s, 0);
    }
    
    return result;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator+(const complex< T >& s)
{
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] + s;
    }
    
    return result;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator-(const complex< T >& s)
{
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] - s;
    }
    
    return result;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator*(const complex< T >& s)
{
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] * s;
    }
    
    return result;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator/(const complex< T >& s)
{
    if (s.re == 0 && s.im == 0)
    {
        uzlmath_error("%s", "division by zero in vector-scalar division.");
    }
    
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] / s;
    }
    
    return result;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator+()
{
    vector< complex< T > > result(size, type);
    memcpy(access::rw(mem), access::rw(result.mem), size * sizeof(complex< T >));
    
    return result;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator-()
{
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = - mem[i];
    }
    
    return result;
}


template< typename T >
inline
vector< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator*(const matrix< T >& mat)
{
    if ((type == vec_type::ROW && size != mat.rows) || (type == vec_type::COLUMN && mat.rows != 1))
    {
        uzlmath_error("%s", "size mismatch in vector-matrix multiplication.");
    }
    
    vector< complex< T > > result(mat.n_cols(), vec_type::ROW);
    
    int M   = (type == vec_type::ROW ? 1 : size);
    int N   = mat.n_cols();
    int K   = (type == vec_type::ROW ? size : 1);
    int LDA = M;
    int LDB = K;
    int LDC = M;
    
    if ( same_type< T, int >::value || same_type< T, short >::value )
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
            result[i].re = static_cast< T >(C[i * 2    ]);
            result[i].im = static_cast< T >(C[i * 2 + 1]);
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
            result[i].re = static_cast< T >(C[i * 2    ]);
            result[i].im = static_cast< T >(C[i * 2 + 1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_mat;
        delete [] C;
    }
    
    return result;
}

template< typename T >
inline
vector< complex< T > > vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator*(const matrix< complex< T > >& mat)
{
    if ((type == vec_type::ROW && size != mat.rows) || (type == vec_type::COLUMN && mat.rows != 1))
    {
        uzlmath_error("%s", "size mismatch in vector-matrix multiplication.");
    }
    
    vector< complex< T > > result(mat.n_cols(), vec_type::ROW);
    
    int M   = (type == vec_type::ROW ? 1 : size);
    int N   = mat.n_cols();
    int K   = (type == vec_type::ROW ? size : 1);
    int LDA = M;
    int LDB = K;
    int LDC = M;
    
    if ( same_type< T, int >::value || same_type< T, short >::value )
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
            result[i].re = static_cast< T >(C[i * 2    ]);
            result[i].im = static_cast< T >(C[i * 2 + 1]);
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
            result[i].re = static_cast< T >(C[i * 2    ]);
            result[i].im = static_cast< T >(C[i * 2 + 1]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_mat;
        delete [] C;
    }
    
    return result;
}

template< typename T >
inline
const vector< complex< T > >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator=(const vector< T >& v)
{
    size = v.size;
    type = v.type;
    
    delete [] mem;
    mem = new complex< T >[size];
    
    if (size > 0)
    {
        size_t i;
        for (i = 0; i < size; ++i)
        {
            mem[i] = complex< T >(v[i], 0);
        }
    }
    
    return *this;
}

template< typename T >
inline
const vector< complex< T > >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator=(const vector< complex< T > >& v)
{
    if ( this == &v )
    {
        return *this;
    }
    
    size = v.size;
    type = v.type;
    
    delete [] mem;
    mem = new complex< T >[size];
    
    if (size > 0)
    {
        memcpy(access::rw(mem), access::rw(v.mem), size * sizeof(complex< T >));
    }
    
    return *this;
}

template< typename T >
inline
const vector< complex< T > >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator=(vector< complex< T > >&& v)
{
    if ( this == &v )
    {
        return *this;
    }
    
    access::rw(size) = v.size;
    access::rw(type) = v.type;
    
    const complex< T >* tmp = mem;
    mem                      = v.mem;
    v.mem                    = tmp;
    
    return *this;
}


template< typename T >
inline
const vector< complex< T > >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator+=(const vector< T >& v)
{
    if (size != size || type != v.type)
    {
        uzlmath_error("%s", "dimension or size mismatch in complex vector-vector multiplication.");
    }
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] += complex< T >(v[i], 0);
    }
    
    return *this;
}

template< typename T >
inline
const vector< complex< T > >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator+=(const vector< complex< T > >& v)
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

template< typename T >
inline
const vector< complex< T > >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator-=(const vector< T >& v)
{
    if (size != size || type != v.type)
    {
        uzlmath_error("%s", "dimension or size mismatch in complex vector-vector multiplication.");
    }
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] -= complex< T >(v[i], 0);
    }
    
    return *this;
}

template< typename T >
inline
const vector< complex< T > >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator-=(const vector< complex< T > >& v)
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

template< typename T >
inline
const vector< complex< T > >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator/=(const vector< T >& v)
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
        
        mem[i] /= complex< T >(v[i], 0);
    }
    
    return *this;
}

template< typename T >
inline
const vector< complex< T > >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator/=(const vector< complex< T > >& v)
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

template< typename T >
inline
const vector< complex< T > >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator%=(const vector< T >& v)
{
    if (type != v.type || size != size)
    {
        uzlmath_error("%s", "type or size mismatch in element-wise complex vector multiplication.");
    }
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] *= complex< T >(v[i], 0);
    }
    
    return *this;
}

template< typename T >
inline
const vector< complex< T > >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator%=(const vector< complex< T > >& v)
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

template< typename T >
inline
const vector< complex< T > >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator*=(const vector< T >& v)
{
    if (type == v.type || type == vec_type::COLUMN || (type == vec_type::ROW && size != v.size))
    {
        uzlmath_error("%s", "size mismatch or wrong vector type in complex vector-vector multiplication.");
    }
    
    size_t i;
    T sum = 0;
    for (i = 0; i < size; ++i)
    {
        sum += mem[i] * complex< T >(v[i], 0);
    }
    
    delete [] mem;
    mem = new complex< T >[1];
    mem[0] = sum;
    
    return *this;
}

template< typename T >
inline
const vector< complex< T > >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator*=(const vector< complex< T > >& v)
{
    if (type == v.type || type == vec_type::COLUMN || (type == vec_type::ROW && size != v.size))
    {
        uzlmath_error("%s", "size mismatch or wrong vector type in complex vector-vector multiplication.");
    }
    
    size_t i;
    T sum = 0;
    for (i = 0; i < size; ++i)
    {
        sum += mem[i] * v[i];
    }
    
    delete [] mem;
    mem = new complex< T >[1];
    mem[0] = sum;
    
    return *this;
}

template< typename T >
inline
const vector< complex< T > >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator+=(const T& s)
{
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] += complex< T >(s, 0);
    }
    
    return *this;
}

template< typename T >
inline
const vector< complex< T > >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator-=(const T& s)
{
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] -= complex< T >(s, 0);
    }
    
    return *this;
}

template< typename T >
inline
const vector< complex< T > >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator*=(const T& s)
{
    
    for (const complex< T >* e = mem; e != mem + size; ++e)
    {
        access::rw(*e) *= complex< T >(s, 0);
    }
    
    return *this;
}

template< typename T >
inline
const vector< complex< T > >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator/=(const T& s)
{
    if (s == 0)
    {
        uzlmath_error("%s", "division by zero in complex vector-scalar division.");
    }
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] /= complex< T >(s, 0);
    }
    
    return *this;
}

template< typename T >
inline
bool vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator==(const vector< T >& v)
{
    bool equal = true;
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        if (mem[i] != complex< T >(v[i], 0))
        {
            equal = false;
            break;
        }
    }
    
    return equal;
}

template< typename T >
inline
bool vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator!=(const vector< T >& v)
{
    return !(*this == v);
}

template< typename T >
inline
bool vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator==(const vector< complex< T > >& v)
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

template< typename T >
inline
bool vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator!=(const vector< complex< T > >& v)
{
    return !(*this == v);
}

template< typename T >
inline
complex< T >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator[](const size_t& idx)
{
    return access::rw(mem[idx]);
}

template< typename T >
inline
constexpr complex< T >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator[](const size_t& idx) const
{
    return access::rw(mem[idx]);
}

template< typename T >
inline
complex< T >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator()(const size_t& idx)
{
    return mem[idx];
}

template< typename T >
inline
constexpr complex< T >& vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::operator()(const size_t& idx) const
{
    return mem[idx];
}

template< typename T >
inline
void vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::ones()
{
    complex< T >* fill_mem = const_cast< complex< T >* >(mem);
    std::fill(fill_mem, fill_mem + size, complex< T >(1, 0));
}

template< typename T >
inline
void vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::zeros()
{
    complex< T >* fill_mem = const_cast< complex< T >* >(mem);
    std::fill(fill_mem, fill_mem + size, complex< T >(0, 0));
}

template< typename T >
inline
void vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::transpose()
{
    type = (type == vec_type::ROW ? vec_type::COLUMN : vec_type::ROW);
}

template< typename T >
inline
void vector< complex< T >, typename if_true< is_num_type< T >::value >::type >::fill(const T& s)
{
    complex< T >* fill_mem = const_cast< complex< T >* >(mem);
    std::fill(fill_mem, fill_mem + size, complex< T >(s, 0));
}



template< typename S >
std::ostream& operator<<(std::ostream& o, const vector< complex< S > >& v)
{
    // setting decimal precesion
    std::ios::fmtflags f( std::cout.flags() );
    o << std::endl;
    
    int width   = 20;
    auto format = std::fixed;
    
    if ( different_type< S, float >::value && different_type< S, double >::value && different_type< S, long double >::value )
    {
        width = 10;
    }
    
    // check values
    size_t i;
    for (i = 0; i < v.size; ++i)
    {
        complex< S > c = v[i];
        
        if (std::abs(c.re) >= 10 || std::abs(c.im) >= 10)
        {
            width   = 22;
            format  = std::fixed;
            
            if ( different_type< S, float >::value && different_type< S, double >::value && different_type< S, long double >::value )
            {
                width = 12;
            }
        }
        
        if (std::abs(c.re) >= 100 || std::abs(c.im) >= 100)
        {
            width   = 24;
            format  = std::fixed;
            
            if ( different_type< S, float >::value && different_type< S, double >::value && different_type< S, long double >::value )
            {
                width = 14;
            }
        }
        
        if (std::abs(c.re) >= 1000 || std::abs(c.im) >= 1000)
        {
            width   = 28;
            format  = std::scientific;
            
            if ( different_type< S, float >::value && different_type< S, double >::value && different_type< S, long double >::value )
            {
                width = 18;
            }
        }
    }
    
    // prepare output and print
    if (v.type == vec_type::ROW)
    {
        for (i = 0; i < v.size; ++i)
        {
            // get entry
            complex< S > c = v[i];
            
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
        o << std::endl;
    }
    else
    {
        for (i = 0; i < v.size; ++i)
        {
            // get entry
            complex< S > c = v[i];
            
            // create string
            std::ostringstream val;
            
            // add real value to string
            val << format << std::setprecision(4) << c.re;
            val << (c.im < 0 ? " - " : " + ") << (c.im == 0 ?  0 : std::abs(c.im)) << "i";
            
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
