//
//  memory_dec.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 29.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_memory_dec_hpp
#define uzlmath_memory_dec_hpp

template< typename eT >
class memory
{
    eT* mem;        //!< Memory that is administrated by this object
    size_t size;    //!< Size of memory allocation
    
public:
    inline memory();
    inline ~memory();
    inline memory(memory&& A);
    inline memory(const memory& A);
    inline memory(const size_t& size);
    
    inline void resize(const size_t& new_size);
    
    inline      size_t n_elements() const;
    
    inline eT& operator[](const size_t& idx);
    inline const eT& operator[](const size_t& idx) const;
    
    inline eT* memptr();
    inline const eT* memptr() const;
};

#endif
