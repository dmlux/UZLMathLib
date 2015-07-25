//
//  array_def.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 29.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_array_def_hpp
#define UZLMathLib_array_def_hpp

UZLMATH_BEGIN

template< typename eT >
inline
array< eT >::array()
    : mem(nullptr)
{}

template< typename eT >
inline
array< eT >::~array()
{
    free(mem);
}

template< typename eT >
inline
array< eT >::array(array&& A)
{
    eT* tmp = mem;
    mem     = A.mem;
    A.mem   = tmp;
    
    size    = A.size;
}

template< typename eT >
inline
array< eT >::array(const array< eT >& A)
{
    memcpy(mem, A.mem, sizeof(eT) * A.size);
    size = A.size;
}

template< typename eT >
inline
array< eT >::array(const size_t& s)
{
    size = s;
    
    if ( void* tmp = malloc(size * sizeof(eT)) )
    {
        mem = static_cast< eT* >(tmp);
    }
    else
    {
        uzlmath_error("%s", "Bad allocation in array allocation.");
    }
}

template< typename eT >
inline
void array< eT >::resize(const size_t& new_size)
{
    size = new_size;
    
    if (new_size == 0)
    {
        free(mem);
        mem = nullptr;
    }
    else
    {
        if (void* tmp = realloc(mem, new_size))
        {
            mem = static_cast< eT* >(tmp);
        }
        else
        {
            uzlmath_error("%s", "Bad allocation in array resize.");
        }
    }
}

template< typename eT >
inline
void array< eT >::fill(const eT& value)
{
    if (is_complex< eT >::value == true)
    {
        std::fill(mem, mem + size, value);
    }
    else
    {
        //memset(mem, value, size * sizeof(eT));
    }
}

template< typename eT >
inline
size_t array< eT >::length() const
{
    return size;
}



template< typename eT >
inline
const eT& array< eT >::operator=(const array< eT >& rhs)
{
    if ( this == &rhs )
    {
        return *this;
    }
    
    size = rhs.size;
    free(mem);
    
    if ( void* tmp = malloc(size * sizeof(eT)) )
    {
        mem = static_cast< eT* >(tmp);
    }
    else
    {
        uzlmath_error("%s", "Bad allocation in array allocation.");
    }
    
    if ( size > 0 )
    {
        memcpy(mem, rhs.mem, size * sizeof(eT));
    }
    
    return *this;
}



template< typename eT >
inline
eT& array< eT >::operator[](const size_t& idx)
{
    return mem[idx];
}

template< typename eT >
inline
const eT& array< eT >::operator[](const size_t& idx) const
{
    return mem[idx];
}

template< typename eT >
inline
eT* array< eT >::memptr()
{
    return mem;
}

template< typename eT >
inline
const eT* array< eT >::memptr() const
{
    return mem;
}

UZLMATH_END

#endif
