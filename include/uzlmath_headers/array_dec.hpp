//
//  array_dec.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 29.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_array_dec_hpp
#define UZLMathLib_array_dec_hpp

UZLMATH_BEGIN

template< typename eT >
class array
{
    eT* mem;        //!< Memory that is administrated by this object
    size_t size;    //!< Size of memory allocation
    
public:
    inline                    array();
    inline                   ~array();
    inline                    array(array&& A);
    inline                    array(const array& A);
    inline                    array(const size_t& s);
    
    inline       void         resize(const size_t& new_size);
    inline       void         fill(const eT& value);
    
    inline       size_t       length() const;
    
    inline const eT&          operator=(const array< eT >& rhs);
    
    inline       eT&          operator[](const size_t& idx);
    inline const eT&          operator[](const size_t& idx) const;
    
    inline       eT*          memptr();
    inline const eT*          memptr() const;
};

UZLMATH_END

#endif
