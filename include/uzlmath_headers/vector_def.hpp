//
//  vector_tpl_def.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 05.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_vector_tpl_def_hpp
#define UZLMathLib_vector_tpl_def_hpp

UZLMATH_BEGIN

/*!
 * @brief           Default constructor for a vector
 * @details         Constructs an empty vector which represents a row vector.
 *                  The memory is not initialized and can contain random
 *                  contents!
 */
template< typename T >
inline
vector< T, if_pod_type< T > >::vector()
    : size(0)
    , type(vec_type::ROW)
    , inj(0)
    , mem(nullptr)
{}

/*!
 * @brief           Constructor for creating a vector of size \f$s\f$
 * @details         Constructs a vector of type \f$t\f$ and size \f$s\f$.
 *                  The memory is not initialized and can contain random
 *                  contents!
 *
 * @param[in]       s The size of the created vector
 * @param[in]       type Type of the created vector where \f$t\f$ can be
 *                  either vec_type::ROW or vec_type::COLUMN
 */
template< typename T >
inline
vector< T, if_pod_type< T > >::vector(const size_t& s, const vec_type& type)
    : inj(0)
    , size(s)
    , type(type)
{
    mem = new T[s];
}

/*!
 * @brief           Constructor for creating a vector of size \f$s\f$ and
 *                  type \f$t\f$ and initial value.
 * @details         Constructs a vector of type \f$t\f$ and size \f$s\f$
 *                  and fills the vector with an given element.
 *
 * @param[in]       s The size of the created vector.
 * @param[in]       initial The initial value for the vector elements
 * @param[in]       type Type of the created vector where \f$t\f$ can be
 *                  either vec_type::ROW or vec_type::COLUMN
 */
template< typename T >
inline
vector< T, if_pod_type< T > >::vector(const size_t& s, const T& initial, const vec_type& type)
    : size(s)
    , type(type)
    , inj(0)
{
    mem = new T[s];
    
    if (s > 0)
    {
        if (initial == 0 || initial == -1)
        {
            memset(mem, initial, size * sizeof(T));
        }
        else
        {
            std::fill(mem, mem + size, initial);
        }
    }
}

/*!
 * @brief           The copy constructor for the vector class.
 * @details         Constructs a new vector that is a copy of the given
 *                  vector. The newly created vector is in the same state
 *                  than the given vector.
 *
 * @param[in]       vec The vector that is supposed to be copied.
 */
template< typename T >
inline
vector< T, if_pod_type< T > >::vector(const vector< T >& vec)
    : size(vec.size)
    , type(vec.type)
    , inj(vec.inj)
{
    mem = new T[size];
    
    if (size > 0)
    {
        memcpy(mem, vec.mem, size * sizeof(T));
    }
}

/*!
 * @brief           The copy constructor for a vector and a new type.
 * @details         Constructs a new vector that is a copy of the given
 *                  vector but with a given new type \f$t\f$. The newly
 *                  created vector is in the same state than the given
 *                  vector.
 *
 * @param[in]       vec The vector that is supposed to be copied.
 * @param[in]       type The new type of the new vector which can either be
 *                  vec_type::ROW or vec_type::COlUMN.
 */
template< typename T >
inline
vector< T, if_pod_type< T > >::vector(const vector< T >& vec, const vec_type& type)
    : size(vec.size)
    , type(type)
    , inj(vec.inj)
{
    mem = new T[size];
    
    if (size > 0)
    {
        memcpy(mem, vec.mem, size * sizeof(T));
    }
}

/*!
 * @brief           The move constructor for a new vector
 * @details         Constructs a new vector by taking the contents of a
 *                  given rvalue vector. The newly created vecotr is in
 *                  the same state than the given vector.
 * 
 * @param[in, out]  vec The vector that is supposed to be copied.
 */
template< typename T >
inline
vector< T, if_pod_type< T > >::vector(vector< T >&& vec)
    : inj(vec.inj)
    , size(vec.size)
    , type(vec.type)
{
    const T* tmp = mem;
    mem           = vec.mem;
    vec.mem       = tmp;
}

/*!
 * @brief           The destructor for the vector.
 * @details         Frees the allocated memory of the vector and preparing
 *                  the object to be released safely.
 */
template< typename T >
inline
vector< T, if_pod_type< T > >::~vector()
{
    delete [] mem;
}



/*!
 * @brief           The addition operator for two vectors.
 * @details         Adds the given vector \f$v\f$ to the current vector and
 *                  creates a new vector object. The addition is done 
 *                  element-wise
 *
 * @param[in]       v The given vector \f$v\f$ which is on the right handside
 *                  of the addition operator.
 *
 * @return          A new vector which contains the element-wise sum of both 
 *                  vectors.
 */
template< typename T >
inline
vector< T > vector< T, if_pod_type< T > >::operator+(const vector< T >& v)
{
    if ( size != v.size || type != v.type)
    {
        uzlmath_error("%s", "size mismatch in vector-vector addition.");
    }
    
    vector< T > result(size, type);
    
    size_t i;
    for (i = 0; i < v.size; ++i)
    {
        result[i] = mem[i] + v[i];
    }
    
    return result;
}

/*!
 * @brief           The subtraction operator for two vectors.
 * @details         Subtracts the given vector \f$v\f$ from the current 
 *                  vector and creates a new vector object. The subtraction
 *                  is done element-wise.
 *
 * @param[in]       v The given vector \f$v\f$ which is on the right handside
 *                  of the addition operator.
 *
 * @return          A new vector which contains the element-wise difference of
 *                  both vectors.
 */
template< typename T >
inline
vector< T > vector< T, if_pod_type< T > >::operator-(const vector< T >& v)
{
    if ( size != v.size || type != v.type)
    {
        uzlmath_error("%s", "size mismatch in vector-vector subtraction.");
    }
    
    vector< T > result(size, type);
    
    size_t i;
    for (i = 0; i < v.size; ++i)
    {
        result[i] = mem[i] + v[i];
    }
    
    return result;
}

/*!
 * @brief           The multiplication operator for two vectors.
 * @details         Multiplies the given vector \f$v\f$ with the current vector and 
 *                  creates a new vector object. The result is either formatted as an 
 *                  \f$1\times 1\f$ or a \f$M\times N\f$ where \f$M,N\in\mathbb{N}^+\f$
 *                  are the sizes of the current vector and the given vector \f$v\f$
 *
 * @param[in]       v The given vector \f$v\f$
 * 
 * @result          The resulting matrix containg the product of both vectors.
 */
template< typename T >
inline
matrix< T > vector< T, if_pod_type< T > >::operator*(const vector< T >& v)
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
    
    matrix< T > result(M, N);
    
    // do multiplication
    if ( same_type< T, int >::value || same_type< T, short >::value )
    {
        
        float* tmp_mem  = new float[size];
        float* tmp_v    = new float[v.size];
        
        size_t i;
        for (i = 0; i < size; ++i)
        {
            tmp_mem[i] = static_cast< float >(mem[i]);
        }
        
        for (i = 0; i < v.size; ++i)
        {
            tmp_v[i]   = static_cast< float >(v[i]);
        }
        
        float* C        = new float[M * N];
        uzl_blas_sgemm(UZLblasNoTrans, UZLblasNoTrans, M, N, K, 1.0, tmp_mem, LDA, tmp_v, LDB, 0.0, C, LDC);
        
        size_t cap_c = M * N;
        for (i = 0; i < cap_c; ++i)
        {
            result.mem[i] = static_cast< T >(C[i]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_v;
        delete [] C;
    }
    else if ( same_type< T, float >::value )
    {
        // Treat pointers as float pointers
        float* A_mem_ptr = reinterpret_cast< float* >(mem);
        float* v_mem_ptr = reinterpret_cast< float* >(v.mem);
        float* C_mem_ptr = reinterpret_cast< float* >(result.mem);
        
        uzl_blas_sgemm(UZLblasNoTrans, UZLblasNoTrans, M, N, K, 1.0, A_mem_ptr, LDA, v_mem_ptr, LDB, 0.0, C_mem_ptr, LDC);
    }
    else if ( same_type< T, double >::value )
    {
        // Treat pointers as double pointers.
        double* A_mem_ptr = reinterpret_cast< double* >(mem);
        double* v_mem_ptr = reinterpret_cast< double* >(v.mem);
        double* C_mem_ptr = reinterpret_cast< double* >(result.mem);
        
        uzl_blas_dgemm(UZLblasNoTrans, UZLblasNoTrans, M, N, K, 1.0, A_mem_ptr, LDA, v_mem_ptr, LDB, 0.0, C_mem_ptr, LDC);
    }
    else
    {
        double* tmp_mem = new double[size];
        double* tmp_v   = new double[v.size];
        
        size_t i;
        for (i = 0; i < size; ++i)
        {
            tmp_mem[i] = static_cast< double >(mem[i]);
        }
        
        for (i = 0; i < v.size; ++i)
        {
            tmp_v[i]   = static_cast< double >(v[i]);
        }
        
        double* C       = new double[M * N];
        uzl_blas_dgemm(UZLblasNoTrans, UZLblasNoTrans, M, N, K, 1.0, tmp_mem, LDA, tmp_v, LDB, 0.0, C, LDC);
        
        size_t cap_c = M * N;
        for (i = 0; i < cap_c; ++i)
        {
            result.mem[i] = static_cast< T >(C[i]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_v;
        delete [] C;
    }
    
    return result;
}

/*!
 * @brief           Element-wise division of two vectors.
 * @details         Computes the element-wise division of two vectors of same type
 *                  and size.
 *
 * @param[in]       v The given vector \f$v\f$
 *
 * @return          The result of element wise division
 */
template< typename T >
inline
vector< T > vector< T, if_pod_type< T > >::operator/(const vector< T >& v)
{
    if (type != v.type || size != size)
    {
        uzlmath_error("%s", "type or size mismatch in element-wise vector division.");
    }
    
    vector< T > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        if (v[i] == 0)
        {
            uzlmath_error("%s", "division by zero in element-wise vector division.");
        }
        
        result[i] = mem[i] / v[i];
    }
    
    return result;
}

/*!
 * @brief           The Hadamard product of two vectors.
 * @details         Computes the element-wise product of two vectors of same type
 *                  and size.
 *
 * @param[in]       v The given vector \f$v\f$
 * 
 * @return          The result of element wise division
 */
template< typename T >
inline
vector< T > vector< T, if_pod_type< T > >::operator%(const vector< T >& v)
{
    if (type != v.type || size != size)
    {
        uzlmath_error("%s", "type or size mismatch in element-wise vector multiplication.");
    }
    
    vector< T > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] * v[i];
    }
    
    return result;
}

/*!
 * @brief           The addition operator for a real and complex vector.
 * @details         Adds the given complex vector \f$v\f$ to the current vector
 *                  and creates a new complex vector. The addition is done
 *                  element-wise.
 * 
 * @param[in]       v The given complex vector \f$v\f$ that is on the right
 *                  handside of the plus operator.
 *
 * @return          A new complex vector containing the result of element-wise
 *                  addition of the current and the given vector.
 */
template< typename T >
inline
vector< complex< T > > vector< T, if_pod_type< T > >::operator+(const vector< complex< T > >& v)
{
    if ( size != v.size || type != v.t)
    {
        uzlmath_error("%s", "size mismatch in vector-vector addition.");
    }
    
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < v.size; ++i)
    {
        result[i] = complex< T >(mem[i], 0) + v[i];
    }
    
    return result;
}

/*!
 * @brief           The subtraction operator for a real and a complex vector.
 * @details         Subtracts the given complex vector \f$v\f$ from the current
 *                  vector and creates a new complex vector object. The 
 *                  subtraction is done element-wise.
 *
 * @param[in]       v The complex given vector \f$v\f$ which is on the right handside
 *                  of the addition operator.
 *
 * @return          A new vector which contains the element-wise difference of
 *                  both vectors.
 */
template< typename T >
inline
vector< complex< T > > vector< T, if_pod_type< T > >::operator-(const vector< complex< T > >& v)
{
    if ( size != v.size || type != v.type)
    {
        uzlmath_error("%s", "size mismatch in vector-vector subtraction.");
    }
    
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < v.size; ++i)
    {
        result[i] = complex< T >(mem[i], 0) - v[i];
    }
    
    return result;
}

/*!
 * @brief           The multiplication operator a real and complex vector.
 * @details         Multiplies the given complex vector \f$v\f$ with the current vector 
 *                  and creates a new complex vector object. The result is either formatted 
 *                  as an \f$1\times 1\f$ or a \f$M\times N\f$ where \f$M,N\in\mathbb{N}^+\f$
 *                  are the sizes of the current vector and the given vector \f$v\f$
 *
 * @param[in]       v The given complex vector \f$v\f$.
 *
 * @result          The resulting matrix containg the product of both vectors.
 */
template< typename T >
inline
matrix< complex< T > > vector< T, if_pod_type< T > >::operator*(const vector< complex< T > >& v)
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

/*!
 * @brief           Element-wise division of real and complex vector.
 * @details         Computes the element-wise division of a real and a complex 
 *                  of same type and size.
 *
 * @param[in]       v The given complex vector \f$v\f$
 *
 * @return          The result of element wise division.
 */
template< typename T >
inline
vector< complex< T > > vector< T, if_pod_type< T > >::operator/(const vector< complex< T > >& v)
{
    if (type != v.type || size != size)
    {
        uzlmath_error("%s", "type or size mismatch in element-wise vector division.");
    }
    
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        if (v[i] == 0)
        {
            uzlmath_error("%s", "division by zero in element-wise vector division.");
        }
        
        result[i] = complex< T >(mem[i], 0) / v[i];
    }
    
    return result;
}

/*!
 * @brief           The Hadamard product of a real and complex vector.
 * @details         Computes the element-wise product of a real and a complex 
 *                  vector of same type and size.
 *
 * @param[in]       v The given vector \f$v\f$
 *
 * @return          The result of element wise division
 */
template< typename T >
inline
vector< complex< T > > vector< T, if_pod_type< T > >::operator%(const vector< complex< T > >& v)
{
    if (type != v.type || size != size)
    {
        uzlmath_error("%s", "type or size mismatch in element-wise vector multiplication.");
    }
    
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = complex< T >(mem[i], 0) * v[i];
    }
    
    return result;
}



/*!
 * @brief           Addition operator for a vector and a scalar.
 * @details         Adds a given scalar value to each element in the vector.
 *
 * @param[in]       s The given scalar value.
 *
 * @return          A new vector containing the sum of the scalar and each 
 *                  vector entry
 */
template< typename T >
inline
vector< T > vector< T, if_pod_type< T > >::operator+(const T& s)
{
    vector< T > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] + s;
    }
    
    return result;
}

/*!
 * @brief           Subtraction operator for a vector and a scalar.
 * @details         Subracts a given scalar value from each element in the
 *                  vector.
 * 
 * @param[in]       s The given scalar value.
 *
 * @return          A new vector containing the differences of each vector element
 *                  and the scalar.
 */
template< typename T >
inline
vector< T > vector< T, if_pod_type< T > >::operator-(const T& s)
{
    vector< T > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] - s;
    }
    
    return result;
}

/*!
 * @brief           Multiplication operator for a vector and a scalar.
 * @details         Multplies a given scalar to each element of the current
 *                  vector.
 *
 * @param[in]       s The given scalar value.
 *
 * @return          A new vector containing the product of each vector element
 *                  and the scalar.
 */
template< typename T >
inline
vector< T > vector< T, if_pod_type< T > >::operator*(const T& s)
{
    vector< T > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] * s;
    }
    
    return result;
}

/*!
 * @brief           Division operator for a vector and a scalar.
 * @details         Divides each element of the current vector by the given
 *                  scalar value.
 *
 * @param[in]       s The given scalar value.
 *
 * @return          A new vector containing the quotient of each vector element
 *                  and the scalar.
 */
template< typename T >
inline
vector< T > vector< T, if_pod_type< T > >::operator/(const T& s)
{
    if (s == 0)
    {
        uzlmath_error("%s", "division by zero in vector-scalar division.");
    }
    
    vector< T > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = mem[i] / s;
    }
    
    return result;
}

/*!
 * @brief           Addition operator for a vector and a complex scalar.
 * @details         Adds the given complex scalar to each element of the
 *                  current vector.
 *
 * @param[in]       s The given complex scalar.
 *
 * @return          A complex vector containing the sum of each element
 *                  and the given scalar.
 */
template< typename T >
inline
vector< complex< T > > vector< T, if_pod_type< T > >::operator+(const complex< T >& s)
{
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = complex< T >(mem[i], 0) + s;
    }
    
    return result;
}

/*!
 * @brief           Subtraction operator for a vector and complex scalar.
 * @details         Subtracts the given complex scalar from each element
 *                  of the current vector.
 *
 * @param[in]       s The given complex scalar.
 *
 * @return          A complex vector containing the difference of each 
 *                  vector element and the scalar value.
 */
template< typename T >
inline
vector< complex< T > > vector< T, if_pod_type< T > >::operator-(const complex< T >& s)
{
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = complex< T >(mem[i], 0) + s;
    }
    
    return result;
}

/*!
 * @brief           Multiplication operator for a vector and a complex scalar.
 * @details         Multiplies the given complex scalar to each element
 *                  in the current vector.
 *
 * @param[in]       s The given complex scalar.
 * 
 * @return          A complex vector containing the product of each
 *                  vector element and the scalar value.
 */
template< typename T >
inline
vector< complex< T > > vector< T, if_pod_type< T > >::operator*(const complex< T >& s)
{
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = complex< T >(mem[i], 0) * s;
    }
    
    return result;
}

/*!
 * @brief           Division operator for a vector and a complex scalar.
 * @details         Divides each element of the current vector by the
 *                  given scalar value.
 *
 * @param[in]       s The given complex scalar value
 *
 * @return
 */
template< typename T >
inline
vector< complex< T > > vector< T, if_pod_type< T > >::operator/(const complex< T >& s)
{
    if (s == 0)
    {
        uzlmath_error("%s", "division by zero in vector-scalar division.");
    }
    
    vector< complex< T > > result(size, type);
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = complex< T >(mem[i], 0) / s;
    }
    
    return result;
}



/*!
 * @brief           Positive sign operator for a vector.
 * @details         returns the vector as it is.
 *
 * @return          returns the vector as it is
 */
template< typename T >
inline
vector< T > vector< T, if_pod_type< T > >::operator+()
{
    vector< T > result;
    memcpy(result.mem, mem, size * sizeof(T));
    
    return result;
}

/*!
 * @brief           Negative sign operator for a vector.
 * @details         Negates each element of the current vector.
 *
 * @return          A vector where each element is negated.
 */
template< typename T >
inline
vector< T > vector< T, if_pod_type< T > >::operator-()
{
    vector< T > result;
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        result[i] = -mem[i];
    }
    
    return result;
}

/*!
 * @brief           Multiplication operator for a vector and a matrix.
 * @details         Multiplies the current vector with the given matrix.
 *
 * @param[in]       mat The given matrix.
 *
 * @return          A vector of type vec_type::ROW containing the product
 *                  of the vector and the matrix.
 */
template< typename T >
inline
vector< T > vector< T, if_pod_type< T > >::operator*(const matrix< T >& mat)
{
    if ((type == vec_type::ROW && size != mat.rows) || (type == vec_type::COLUMN && mat.rows != 1))
    {
        uzlmath_error("%s", "size mismatch in vector-matrix multiplication.");
    }
    
    vector< T > result(mat.n_cols(), vec_type::ROW);
    
    int M   = (type == vec_type::ROW ? 1 : size);
    int N   = mat.n_cols();
    int K   = (type == vec_type::ROW ? size : 1);
    int LDA = M;
    int LDB = K;
    int LDC = M;
    
    if ( same_type< T, int >::value || same_type< T, short >::value )
    {
        float* tmp_mem  = new float[size];
        float* tmp_mat  = new float[mat.n_cols() * mat.n_rows()];
        
        size_t i;
        for (i = 0; i < size; ++i)
        {
            tmp_mem[i] = static_cast< float >(mem[i]);
        }
        
        for (i = 0; i < mat.n_cols() * mat.n_rows(); ++i)
        {
            tmp_mat[i] = static_cast< float >(mat.mem[i]);
        }
        
        float* C        = new float[M * N];
        uzl_blas_sgemm(UZLblasNoTrans, UZLblasNoTrans, M, N, K, 1.0, tmp_mem, LDA, tmp_mem, LDB, 0.0, C, LDC);
        
        size_t cap_c = M * N;
        for (i = 0; i < cap_c; ++i)
        {
            result[i] = static_cast< T >(C[i]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_mat;
        delete [] C;
    }
    else if ( same_type< T, float >::value )
    {
        float* A_mem_ptr = reinterpret_cast< float* >(mem);
        float* B_mem_ptr = reinterpret_cast< float* >(mat.mem);
        float* C_mem_ptr = reinterpret_cast< float* >(result.mem);
        
        uzl_blas_sgemm(UZLblasNoTrans, UZLblasNoTrans, M, N, K, 1.0, A_mem_ptr, LDA, B_mem_ptr, LDB, 0.0, C_mem_ptr, LDC);
    }
    else if ( same_type< T, double >::value )
    {
        double* A_mem_ptr = reinterpret_cast< double* >(mem);
        double* B_mem_ptr = reinterpret_cast< double* >(mat.mem);
        double* C_mem_ptr = reinterpret_cast< double* >(result.mem);
        
        uzl_blas_dgemm(UZLblasNoTrans, UZLblasNoTrans, M, N, K, 1.0, A_mem_ptr, LDA, B_mem_ptr, LDB, 0.0, C_mem_ptr, LDC);
    }
    else
    {
        double* tmp_mem = new double[size];
        double* tmp_mat = new double[mat.n_cols() * mat.n_rows()];
        
        size_t i;
        for (i = 0; i < size; ++i)
        {
            tmp_mem[i] = static_cast< double >(mem[i]);
        }
        
        for (i = 0; i < mat.n_cols() * mat.n_rows(); ++i)
        {
            tmp_mat[i] = static_cast< double >(mat.mem[i]);
        }
        
        double* C       = new double[M * N];
        uzl_blas_dgemm(UZLblasNoTrans, UZLblasNoTrans, M, N, K, 1.0, tmp_mem, LDA, tmp_mem, LDB, 0.0, C, LDC);
        
        size_t cap_c = M * N;
        for (i = 0; i < cap_c; ++i)
        {
            result[i] = static_cast< T >(C[i]);
        }
        
        delete [] tmp_mem;
        delete [] tmp_mat;
        delete [] C;
    }
    
    return result;
}



/*!
 * @brief           The greater than operator for two vectors.
 * @details         Checks if each element of the current vector is
 *                  greater than the corresponding element in the 
 *                  given vector.
 *                  \f[
 *                      \left[\overrightarrow{v}\right]_k > \left[\overrightarrow{w}\right]_k
 *                  \f]
 *
 * @param[in]       v The vector that is supposed to be checked
 *
 * @return          True if each vector element is greater than the
 *                  corresponding element in the given vector, false
 *                  else.
 */
template< typename T >
inline
bool vector< T, if_pod_type< T > >::operator>(const vector< T >& v)
{
    if (size != size)
    {
        uzlmath_error("%s", "size mismatch in greater than vector comparison.");
    }
    
    bool greater = true;
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        if (mem[i] <= v[i])
        {
            greater = false;
            break;
        }
    }
    
    return greater;
}

/*!
 * @brief           The lower than operator for two vectors.
 * @details         Checks if each element of the current vector is
 *                  lower than the corresponding element in the
 *                  given vector.
 *                  \f[
 *                      \left[\overrightarrow{v}\right]_k < \left[\overrightarrow{w}\right]_k
 *                  \f]
 *
 * @param[in]       v The vector that is supposed to be checked
 *
 * @return          True if each vector element is lower than the
 *                  corresponding element in the given vector, false
 *                  else.
 */
template< typename T >
inline
bool vector< T, if_pod_type< T > >::operator<(const vector< T >& v)
{
    if (size != size)
    {
        uzlmath_error("%s", "size mismatch in lower than vector comparison.");
    }
    
    bool lower = true;
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        if (mem[i] >= v[i])
        {
            lower = false;
            break;
        }
    }
    
    return lower;
}



/*!
 * @brief           The assignment operator for a vector.
 * @details         Assigns the given vector to the current vector.
 *                  The newly created vector is in the same state than
 *                  the given vector.
 *
 * @param[in]       v The given vector that is supposed to be assignend.
 *
 * @return          A reference to the current vector.
 */
template< typename T >
inline
const vector< T >& vector< T, if_pod_type< T > >::operator=(const vector< T >& v)
{
    if ( this == &v )
    {
        return *this;
    }
    
    size = v.size;
    type = v.type;
    
    delete [] mem;
    mem = new T[size];
    
    if (size > 0)
    {
        memcpy(mem, v.mem, size * sizeof(T));
    }
        
    return *this;
}

/*!
 * @brief           The assigment move operator.
 * @details         Assigns the contents of a given vector \f$v\f$
 *                  which is a rvalue vector, by moving its contents
 *                  to the current vector object.
 */
template< typename T >
inline
const vector< T >& vector< T, if_pod_type< T > >::operator=(vector< T >&& v)
{
    if ( this == &v )
    {
        return *this;
    }
    
    size = v.size;
    type = v.type;
    
    T* tmp = mem;
    mem     = v.mem;
    v.mem   = tmp;
    
    return *this;
}



/*!
 * @brief           The addition assignment operator for two vectors.
 * @details         Adds the given vector \f$v\f$ to the current vector
 *                  and returning a reference to the current vector.
 *
 * @param[in]       v The given vector.
 *
 * @return          A reference to the current vector object.
 */
template< typename T >
inline
const vector< T >& vector< T, if_pod_type< T > >::operator+=(const vector< T >& v)
{
    if (size != size || type != v.type)
    {
        uzlmath_error("%s", "dimension or size mismatch in vector-vector multiplication.");
    }
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] += v[i];
    }
    
    return *this;
}

/*!
 * @brief           The subtraction assignment operator for two vectors.
 * @details         Subtracts the given vector \f$v\f$ from the current
 *                  vector by subtracting the elements element-wise.
 *
 * @param[in]       v The given vector.
 *
 * @return          A reference to the current vector object.
 */
template< typename T >
inline
const vector< T >& vector< T, if_pod_type< T > >::operator-=(const vector< T >& v)
{
    if (size != v.size || type != v.type)
    {
        uzlmath_error("%s", "dimension or size mismatch in vector-vector subtraction.");
    }
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] -= v[i];
    }
    
    return *this;
}

/*!
 * @brief           The element-wise division operator for two vectors.
 * @details         Divides the current vector by the given vector element-
 *                  wise.
 *
 * @param[in]       v The given vector
 *
 * @return          A reference to the current vector object.
 */
template< typename T >
inline
const vector< T >& vector< T, if_pod_type< T > >::operator/=(const vector< T >& v)
{
    if (type != v.type || size != size)
    {
        uzlmath_error("%s", "type or size mismatch in element-wise vector division.");
    }
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        if (v[i] == 0)
        {
            uzlmath_error("%s", "division by zero in element-wise vector division.");
        }
        
        mem[i] /= v[i];
    }
    
    return *this;
}

/*!
 * @brief           The Hadamard product for two vectors.
 * @details         Performing an element-wise multiplication of the current
 *                  and the given vector.
 *
 * @param[in]       v The given vector
 *
 * @return          A reference to the current vector object.
 */
template< typename T >
inline
const vector< T >& vector< T, if_pod_type< T > >::operator%=(const vector< T >& v)
{
    if (type != v.type || size != size)
    {
        uzlmath_error("%s", "type or size mismatch in element-wise vector multiplication.");
    }
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] *= v[i];
    }
    
    return *this;
}



/*!
 * @brief           The addition assignment operator for a vector and a scalar.
 * @details         Adds a given scalar to the current vector. The scalar value
 *                  gets added to each vector element.
 *
 * @param[in]       s The scalar that is supposed to be added to the current vector.
 * 
 * @return          The reference to the current vector.
 */
template< typename T >
inline
const vector< T >& vector< T, if_pod_type< T > >::operator+=(const T& s)
{
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] += s;
    }
    
    return *this;
}

/*!
 * @brief           The subtraction assignment operator for a vector and a scalar.
 * @details         Subtracts a given scalar from the current vector. The scalar
 *                  value gets subtracted from each vector element.
 *
 * @param[in]       s The scalar that is supposed to be subtracted from the current
 *                  vector.
 *
 * @return          The reference to the current vector.
 */
template< typename T >
inline
const vector< T >& vector< T, if_pod_type< T > >::operator-=(const T& s)
{
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] -= s;
    }
    
    return *this;
}

/*!
 * @brief           The multiplication assignment operator for a vector and a scalar.
 * @details         Multiplies a given scalar to the current vector. The scalar value
 *                  gets multiplied to each vector element.
 *
 * @param[in]       s The scalar that is supposed to be multiplied to the current vector.
 *
 * @return          The reference to the current vector.
 */
template< typename T >
inline
const vector< T >& vector< T, if_pod_type< T > >::operator*=(const T& s)
{
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] *= s;
    }
    
    return *this;
}

/*!
 * @brief           The division assignment operator for a vector and a scalar.
 * @details         Divides a given scalar from the current vector. Each element of the
 *                  vector gets divided by the given scalar.
 *
 * @param[in]       s The scalar by which each element of the current vector gets divided.
 * 
 * @return          The reference to the current vector.
 */
template< typename T >
inline
const vector< T >& vector< T, if_pod_type< T > >::operator/=(const T& s)
{
    if (s == 0)
    {
        uzlmath_error("%s", "division by zero in vector-scalar division.");
    }
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        mem[i] /= s;
    }
    
    return *this;
}



/*!
 * @brief           The equal comparison operator for two vectors.
 * @details         Compares the two vectors for equality.
 *
 * @param[in]       v The vector that is supposed to be compared to the current vector.
 *
 * @return          True if both vectors are element-wise equal, else false.
 */
template< typename T >
inline
bool vector< T, if_pod_type< T > >::operator==(const vector< T >& v)
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

/*!
 * @brief           The unequal comparison operator for two vectors.
 * @details         Compares the two vectors for unequality.
 *
 * @param[in]       v The vector that is supposed to be compared to the current vector.
 *
 * @return          True if both vectors are element-wise unequal, else true.
 */
template< typename T >
inline
bool vector< T, if_pod_type< T > >::operator!=(const vector< T >& v)
{
    return !(*this == v);
}

/*!
 * @brief           The greater equals operator for two vectors.
 * @details         Checks if each element of the current vector is
 *                  greater or equal to the corresponding element in the
 *                  given vector.
 *                  \f[
 *                      \left[\overrightarrow{v}\right]_k \geq \left[\overrightarrow{w}\right]_k
 *                  \f]
 *
 * @param[in]       v The vector that is supposed to be checked
 *
 * @return          True if each vector element is greater or equal the
 *                  corresponding element in the given vector, false
 *                  else.
 */
template< typename T >
inline
bool vector< T, if_pod_type< T > >::operator>=(const vector< T >& v)
{
    if (size != size)
    {
        uzlmath_error("%s", "size mismatch in greater equals vector comparison.");
    }
    
    bool geq = true;
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        if (mem[i] < v[i])
        {
            geq = false;
            break;
        }
    }
    
    return geq;
}

/*!
 * @brief           The lower equuals operator for two vectors.
 * @details         Checks if each element of the current vector is
 *                  lower or equal to the corresponding element in the
 *                  given vector.
 *                  \f[
 *                      \left[\overrightarrow{v}\right]_k \leq \left[\overrightarrow{w}\right]_k
 *                  \f]
 *
 * @param[in]       v The vector that is supposed to be checked
 *
 * @return          True if each vector element is lower or equal to the
 *                  corresponding element in the given vector, false
 *                  else.
 */
template< typename T >
inline
bool vector< T, if_pod_type< T > >::operator<=(const vector< T >& v)
{
    if (size != size)
    {
        uzlmath_error("%s", "size mismatch in lower than vector comparison.");
    }
    
    bool leq = true;
    
    size_t i;
    for (i = 0; i < size; ++i)
    {
        if (mem[i] > v[i])
        {
            leq = false;
            break;
        }
    }
    
    return leq;
}



/*!
 * @brief           The write subscript operator for the vector object.
 * @details         Overload for the subscript operator. The subscript operator is used 
 *                  to access the elements of the vector by using square brackets on the 
 *                  vector object like this.
 *                  @code
 *                      vector<double> v(5, vec_type::ROW);
 *                      v.fill(5);
 *
 *                      // use subscript to overwrite the third value
 *                      v[3] = 7;
 *                      
 *                      // The content of v looks like this
 *                      // v =
 *                      //     5  5  7  5  5
 *                  @endcode
 *                  
 * @param[in]       idx The index of the element that is supposed to be returned.
 * 
 * @return          The reference to element at index idx
 */
template< typename T >
inline
T& vector< T, if_pod_type< T > >::operator[](const size_t& idx)
{
    return access::rw(mem[idx]);
}

/*!
 * @brief           The read subscript operator for the vector objects.
 * @details         Overload for the subscript operator. The subscript operator is used
 *                  to access the elements of the vector by using square brackets on the
 *                  vector object like this.
 *                  @code
 *                      double a = 3;
 *                      vector<double> v(5, vec_type::ROW)
 *                      v.fill(5);
 *
 *                      // use subscript to read from the vector
 *                      a = v[5];
 *                  
 *                      // The content of a looks like this
 *                      // a = 5
 *                  @endcode
 *
 * @param[in]       idx The index of the element that is supposed to be read.
 *
 * @return          The element at index idx
 */
template< typename T >
inline
const T& vector< T, if_pod_type< T > >::operator[](const size_t& idx) const
{
    return access::rw(mem[idx]);
}



/*!
 * @brief           The write access operator for the vector object.
 * @details         Overload for the access operator. The access operator is used
 *                  to access the elements of the vector by using parenthesis on the
 *                  vector object like this.
 *                  @code
 *                      vector<double> v(5, vec_type::ROW);
 *                      v.fill(5);
 *
 *                      // use subscript to overwrite the third value
 *                      v(3) = 7;
 *
 *                      // The content of v looks like this
 *                      // v =
 *                      //     5  5  7  5  5
 *                  @endcode
 *
 * @param[in]       idx The index of the element that is supposed to be returned.
 *
 * @return          The reference to element at index idx
 */
template< typename T >
inline
T& vector< T, if_pod_type< T > >::operator()(const size_t& idx)
{
    return mem[idx];
}

/*!
 * @brief           The read access operator for the vector objects.
 * @details         Overload for the access operator. The access operator is used
 *                  to access the elements of the vector by using parenthesis on the
 *                  vector object like this.
 *                  @code
 *                      double a = 3;
 *                      vector<double> v(5, vec_type::ROW)
 *                      v.fill(5);
 *
 *                      // use subscript to read from the vector
 *                      a = v(5);
 *
 *                      // The content of a looks like this
 *                      // a = 5
 *                  @endcode
 *
 * @param[in]       idx The index of the element that is supposed to be read.
 *
 * @return          The element at index idx
 */
template< typename T >
inline
const T& vector< T, if_pod_type< T > >::operator()(const size_t& idx) const
{
    return mem[idx];
}



/*!
 * @brief           Filling the vector with one on each entry.
 * @details         Each entry of the current vector gets overwritten with 1.
 */
template< typename T >
inline
void vector< T, if_pod_type< T > >::ones()
{
    std::fill(mem, mem + size, 1);
}

/*!
 * @brief           Filling the vector with zeros on each entry.
 * @details         Each entry of the current vector gets overwritten with 0.
 */
template< typename T >
inline
void vector< T, if_pod_type< T > >::zeros()
{
    memset(mem, 0, size * sizeof(T));
}

/*!
 * @brief           Transposes the current vector.
 * @details         Switching vector type from vec_type::ROW to vec_type::COLUMN or
 *                  from vec_type::COLUMN to vec_type::ROW
 */
template< typename T >
inline
void vector< T, if_pod_type< T > >::transpose()
{
    if (type == vec_type::ROW)
    {
        type = vec_type::COLUMN;
    }
    else
    {
        type = vec_type::ROW;
    }
}

/*!
 * @brief           Filling the vector with given scalar value.
 * @details         Each entry of the current vector gets overwritten with the
 *                  given scalar value.
 */
template< typename T >
inline
void vector< T, if_pod_type< T > >::fill(const T& s)
{
    if (s == 0 || s == -1)
    {
        memset(mem, s, size * sizeof(T));
    }
    else
    {
        std::fill(mem, mem + size, s);
    }
}



/*!
 * @brief           Outstream operator overload.
 * @details         The out-steam operator is used to print the vector
 *                  in a nice form over the std::cout stream.
 *
 * @param[in,out]   o The stream object
 * @param[in]       v The vector that should be printed
 *
 * @return          The reference to the given out-stream.
 */
template< typename S >
std::ostream& operator<<(std::ostream& o, const vector< S >& v)
{
    // setting decimal precesion
    std::ios::fmtflags f( std::cout.flags() );
    o << std::endl;
    
    int width   = 10;
    auto format = std::fixed;
    
    // reduce size for integers
    if ( different_type< S, float >::value && different_type< S, double >::value && different_type< S, long double >::value )
    {
        width = 5;
    }
    
    // check values
    size_t i;
    for (i = 0; i < v.size; ++i)
    {
        S val = v[i];
        if (UZL_ABS(val) >= 10)
        {
            width   = 11;
            format  = std::fixed;
            
            if ( different_type< S, float >::value && different_type< S, double >::value && different_type< S, long double >::value )
            {
                width = 6;
            }
        }
        
        if (UZL_ABS(val) >= 100)
        {
            width   = 12;
            format  = std::fixed;
            
            if ( different_type< S, float >::value && different_type< S, double >::value && different_type< S, long double >::value )
            {
                width = 7;
            }
        }
        
        if (UZL_ABS(val) >= 1000)
        {
            width   = 14;
            format  = std::scientific;
            
            if ( different_type< S, float >::value && different_type< S, double >::value && different_type< S, long double >::value )
            {
                width = 10;
            }
        }
    }
    
    if (v.vecType() == vec_type::ROW)
    {
        for (i = 0; i < v.size; ++i)
        {
            // get entry
            S val = v[i];
            
            // create string
            o << std::setfill(' ');
            o << std::right << std::setw(width);
            o << format << std::setprecision(4) << val;
        }
        o << std::endl;
    }
    else
    {
        for (i = 0; i < v.size; ++i)
        {
            // get entry
            S val = v[i];
            
            // create string
            o << std::setfill(' ');
            o << std::right << std::setw(width);
            o << format << std::setprecision(4) << val << std::endl;
        }
    }
    
    std::cout.flags( f );
    return o;
}

UZLMATH_END

#endif /* vector_def.hpp */
