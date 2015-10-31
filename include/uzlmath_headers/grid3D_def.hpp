//
//  grid3D_def.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 19.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_grid3D_def_hpp
#define UZLMathLib_grid3D_def_hpp

UZLMATH_BEGIN

template< typename T >
inline
grid3D< T, typename if_true< is_num_type< T >::value >::type >::grid3D()
    : rows(0)
    , cols(0)
    , lays(0)
    , mem(nullptr)
{}

template< typename T >
inline
grid3D< T, typename if_true< is_num_type< T >::value >::type >::grid3D(const size_t& rows, const size_t& cols, const size_t& lays)
    : rows(rows)
    , cols(cols)
    , lays(lays)
{
    mem = new T[rows * cols * lays];
}

template< typename T >
inline
grid3D< T, typename if_true< is_num_type< T >::value >::type >::grid3D(const size_t& rcl)
    : rows(rcl)
    , cols(rcl)
    , lays(rcl)
{
    mem = new T[rows * cols * lays];
}

template< typename T >
inline
grid3D< T, typename if_true< is_num_type< T >::value >::type >::grid3D(const size_t& rows, const size_t& cols, const size_t& lays, const T& initial)
    : rows(rows)
    , cols(cols)
    , lays(lays)
{
    size_t cap  = rows * cols * lays;
    mem         = new T[cap];
    
    if (cap > 0)
    {
        if (initial == 0 || initial == -1)
        {
            memset(mem, initial, cap * sizeof(T));
        }
        else
        {
            std::fill(mem, mem + cap, initial);
        }
    }
}

template< typename T >
inline
grid3D< T, typename if_true< is_num_type< T >::value >::type >::grid3D(const size_t& rcl, const T& initial)
    : rows(rcl)
    , cols(rcl)
    , lays(rcl)
{
    size_t cap  = rows * cols * lays;
    mem         = new T[cap];
    
    if (cap > 0)
    {
        if (initial == 0 || initial == -1)
        {
            memset(mem, initial, cap * sizeof(T));
        }
        else
        {
            std::fill(mem, mem + cap, initial);
        }
    }
}

template< typename T >
inline
grid3D< T, typename if_true< is_num_type< T >::value >::type >::grid3D(const grid3D< T >& c)
    : rows(c.rows)
    , cols(c.cols)
    , lays(c.lays)
{
    size_t cap  = rows * cols * lays;
    mem         = new T[cap];
    
    if (cap > 0)
    {
        memcpy(mem, c.mem, c.rows * c.cols * c.lays * sizeof(T));
    }
}

template< typename T >
inline
grid3D< T, typename if_true< is_num_type< T >::value >::type >::grid3D(grid3D< T >&& c)
    : rows(c.rows)
    , cols(c.cols)
    , lays(c.lays)
{
    T* tmp = mem;
    mem     = c.mem;
    c.mem   = tmp;
}

template< typename T >
inline
grid3D< T, typename if_true< is_num_type< T >::value >::type >::~grid3D()
{
    delete [] mem;
}

template< typename T >
inline
const grid3D< T >& grid3D< T, typename if_true< is_num_type< T >::value >::type >::operator=(const grid3D< T >& c)
{
    if ( *this == &c )
    {
        return *this;
    }
    
    rows = c.rows;
    cols = c.cols;
    lays = c.lays;
    
    delete [] mem;
    
    size_t cap  = rows * cols * lays;
    mem         = new T[cap];
    
    if (cap > 0)
    {
        memcpy(mem, c.mem, cap * sizeof(T));
    }
    
    return *this;
}

template< typename T >
inline
const grid3D< T >& grid3D< T, typename if_true< is_num_type< T >::value >::type >::operator=(grid3D< T >&& c)
{
    if ( *this == &c )
    {
        return *this;
    }
    
    rows = c.rows;
    cols = c.cols;
    lays = c.lays;
    
    T* tmp = mem;
    mem     = c.mem;
    c.mem   = tmp;
    
    return *this;
}

template< typename T >
inline
T& grid3D< T, typename if_true< is_num_type< T >::value >::type >::operator()(const size_t& r, const size_t& c, const size_t& l)
{
    return mem[rows * cols * l + cols * c + r];
}

template< typename T >
inline
const T& grid3D< T, typename if_true< is_num_type< T >::value >::type >::operator()(const size_t& r, const size_t& c, const size_t& l) const
{
    return mem[rows * cols * l + cols * c + r];
}



template< typename S >
std::ostream& operator<<(std::ostream& o, const grid3D< S >& c)
{
    std::ios::fmtflags f( std::cout.flags() );
    o << std::endl;
    
    int width   = 10;
    auto format = std::fixed;
    
    // reduce size for integers
    if ( different_type< S, float >::value && different_type< S, double >::value && different_type< S, long double >::value )
    {
        width = 5;
    }
    
    // check values for nice formatting
    size_t x, y, z;
    for (z = 0; z < c.lays; ++z)
    {
        for (x = 0; x < c.rows; ++x)
        {
            for (y = 0; y < c.cols; ++y)
            {
                S val = c(x, y, z);
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
        }
    }
    
    // setting decimal precesion
    for (z = 0; z < c.lays; ++z)
    {
        o << "layer[" << z << "]" << std::endl;
        for (x = 0; x < c.rows; ++x)
        {
            for (y = 0; y < c.cols; ++y)
            {
                // get entry
                S val = c(x, y, z);
                
                // create string
                o << std::setfill(' ');
                o << std::right << std::setw(width);
                o << format << std::setprecision(4) << val;
            }
            o << std::endl;
        }
        o << std::endl;
    }
    
    std::cout.flags( f );
    return o;
}

UZLMATH_END

#endif /* grid3D_def.hpp */
