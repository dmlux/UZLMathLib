//
//  memory_def.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 29.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_memory_def_hpp
#define uzlmath_memory_def_hpp

template< typename eT >
inline
memory< eT >::memory()
    : mem(nullptr)
{}

template< typename eT >
inline
memory< eT >::~memory()
{
    free(mem);
}

template< typename eT >
inline
memory< eT >::memory(memory&& A)
{
    eT* tmp = mem;
    mem     = A.mem;
    A.mem   = tmp;
    
    size = A.size;
}

template< typename eT >
inline
memory< eT >::memory(const memory< eT >& A)
{
    memcpy(mem, A.mem, sizeof(eT) * A.size);
    size = A.size;
}

template< typename eT >
inline
memory< eT >::memory(const size_t& size)
{
    this->size = size;
    
    if (void* tmp = malloc(size * sizeof(eT)))
    {
        mem = static_cast< eT* >(tmp);
    }
    else
    {
        std::cout << "** uzlmath error: Bad allocation in memory allocation. **" << std::endl;
        exit(EXIT_FAILURE);
    }
}

template< typename eT >
inline
void memory< eT >::resize(const size_t& new_size)
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
            std::cout << "** uzlmath error: Bad allocation in memory resize. **" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}

template< typename eT >
inline
size_t memory< eT >::n_elements() const
{
    return size;
}

template< typename eT >
inline
eT& memory< eT >::operator[](const size_t& idx)
{
    return mem[idx];
}

template< typename eT >
inline
const eT& memory< eT >::operator[](const size_t& idx) const
{
    return mem[idx];
}

template< typename eT >
inline
eT* memory< eT >::memptr()
{
    return mem;
}

template< typename eT >
inline
const eT* memory< eT >::memptr() const
{
    return mem;
}

#endif
