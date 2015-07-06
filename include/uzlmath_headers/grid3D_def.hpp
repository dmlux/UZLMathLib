//
//  grid3D_def.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 19.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef uzlmath_grid3D_def_hpp
#define uzlmath_grid3D_def_hpp

UZLMATH_BEGIN

template< typename eT >
inline
grid3D< eT >::grid3D()
    : rows(0)
    , cols(0)
    , lays(0)
    , mem(nullptr)
{}

template< typename eT >
inline
grid3D< eT >::grid3D(const size_t& rows, const size_t& cols, const size_t& lays)
    : rows(rows)
    , cols(cols)
    , lays(lays)
{
    mem = new eT[rows * cols * lays];
}

template< typename eT >
inline
grid3D< eT >::grid3D(const size_t& rcl)
    : rows(rcl)
    , cols(rcl)
    , lays(rcl)
{
    mem = new eT[rows * cols * lays];
}

template< typename eT >
inline
grid3D< eT >::grid3D(const size_t& rows, const size_t& cols, const size_t& lays, const eT& initial)
    : rows(rows)
    , cols(cols)
    , lays(lays)
{
    size_t cap  = rows * cols * lays;
    mem         = new eT[cap];
    
    if (cap > 0)
    {
        if (initial == 0 || initial == -1)
        {
            memset(mem, initial, cap * sizeof(eT));
        }
        else
        {
            std::fill(mem, mem + cap, initial);
        }
    }
}

template< typename eT >
inline
grid3D< eT >::grid3D(const size_t& rcl, const eT& initial)
    : rows(rcl)
    , cols(rcl)
    , lays(rcl)
{
    size_t cap  = rows * cols * lays;
    mem         = new eT[cap];
    
    if (cap > 0)
    {
        if (initial == 0 || initial == -1)
        {
            memset(mem, initial, cap * sizeof(eT));
        }
        else
        {
            std::fill(mem, mem + cap, initial);
        }
    }
}

template< typename eT >
inline
grid3D< eT >::grid3D(const grid3D< eT >& c)
    : rows(c.rows)
    , cols(c.cols)
    , lays(c.lays)
{
    size_t cap  = rows * cols * lays;
    mem         = new eT[cap];
    
    if (cap > 0)
    {
        memcpy(mem, c.mem, c.rows * c.cols * c.lays * sizeof(eT));
    }
}

template< typename eT >
inline
grid3D< eT >::grid3D(grid3D< eT >&& c)
    : rows(c.rows)
    , cols(c.cols)
    , lays(c.lays)
{
    eT* tmp = mem;
    mem     = c.mem;
    c.mem   = tmp;
}

template< typename eT >
inline
grid3D< eT >::~grid3D()
{
    delete [] mem;
}

template< typename eT >
inline
const grid3D< eT >& grid3D< eT >::operator=(const grid3D< eT >& c)
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
    mem         = new eT[cap];
    
    if (cap > 0)
    {
        memcpy(mem, c.mem, cap * sizeof(eT));
    }
    
    return *this;
}

template< typename eT >
inline
const grid3D< eT >& grid3D< eT >::operator=(grid3D< eT >&& c)
{
    if ( *this == &c )
    {
        return *this;
    }
    
    rows = c.rows;
    cols = c.cols;
    lays = c.lays;
    
    eT* tmp = mem;
    mem     = c.mem;
    c.mem   = tmp;
    
    return *this;
}

template< typename eT >
inline
eT& grid3D< eT >::operator()(const size_t& r, const size_t& c, const size_t& l)
{
    return mem[rows * cols * l + cols * c + r];
}

template< typename eT >
inline
const eT& grid3D< eT >::operator()(const size_t& r, const size_t& c, const size_t& l) const
{
    return mem[rows * cols * l + cols * c + r];
}

template< typename eT >
inline
size_t grid3D< eT >::n_rows() const
{
    return rows;
}

template< typename eT >
inline
size_t grid3D< eT >::n_cols() const
{
    return cols;
}

template< typename eT >
inline
size_t grid3D< eT >::n_lays() const
{
    return lays;
}

template< typename eT >
inline
eT* grid3D< eT >::memptr()
{
    return mem;
}

template< typename eT >
inline
const eT* grid3D< eT >::memptr() const
{
    return mem;
}



template< typename S >
std::ostream& operator<<(std::ostream& o, const grid3D< S >& c)
{
    std::ios::fmtflags f( std::cout.flags() );
    o << std::endl;
    
    int width   = 10;
    auto format = std::fixed;
    
    // reduce size for integers
    if (is_float< S >::value == false && is_double< S >::value == false && is_ldouble< S >::value == false)
    {
        width = 5;
    }
    
    // check values for nice formatting
    size_t x, y, z;
    for (z = 0; z < c.n_lays(); ++z)
    {
        for (x = 0; x < c.n_rows(); ++x)
        {
            for (y = 0; y < c.n_cols(); ++y)
            {
                S val = c(x, y, z);
                if (UZL_ABS(val) >= 10)
                {
                    width   = 11;
                    format  = std::fixed;
                    
                    if (is_float< S >::value == false && is_double< S >::value == false && is_ldouble< S >::value == false)
                    {
                        width = 6;
                    }
                }
                
                if (UZL_ABS(val) >= 100)
                {
                    width   = 12;
                    format  = std::fixed;
                    
                    if (is_float< S >::value == false && is_double< S >::value == false && is_ldouble< S >::value == false)
                    {
                        width = 7;
                    }
                }
                
                if (UZL_ABS(val) >= 1000)
                {
                    width   = 14;
                    format  = std::scientific;
                    
                    if (is_float< S >::value == false && is_double< S >::value == false && is_ldouble< S >::value == false)
                    {
                        width = 10;
                    }
                }
            }
        }
    }
    
    // setting decimal precesion
    for (z = 0; z < c.n_lays(); ++z)
    {
        o << "layer[" << z << "]" << std::endl;
        for (x = 0; x < c.n_rows(); ++x)
        {
            for (y = 0; y < c.n_cols(); ++y)
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

#endif
