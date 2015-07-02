//
//  cx_grid3D_def.hpp
//  uzlmath
//
//  Created by Denis-Michael Lux on 08.06.15.
//
//

#ifndef uzlmath_cx_grid3D_def_h
#define uzlmath_cx_grid3D_def_h

template< typename eT >
inline
grid3D< complex< eT > >::grid3D()
    : rows(0)
    , cols(0)
    , lays(0)
    , mem(nullptr)
{}

template< typename eT >
inline
grid3D< complex< eT > >::grid3D(const size_t& rows, const size_t& cols, const size_t& lays)
    : rows(rows)
    , cols(cols)
    , lays(lays)
{
    mem = new complex< eT >[rows * cols * lays];
}

template< typename eT >
inline
grid3D< complex< eT > >::grid3D(const size_t& rcl)
    : rows(rcl)
    , cols(rcl)
    , lays(rcl)
{
    mem = new complex< eT >[rows * cols * lays];
}

template< typename eT >
inline
grid3D< complex< eT > >::grid3D(const size_t& rows, const size_t& cols, const size_t& lays, const complex< eT >& initial)
    : rows(rows)
    , cols(cols)
    , lays(lays)
{
    mem = new complex< eT >[rows * cols * lays];
    
    std::fill(mem, mem + rows * cols * lays, initial);
}

template< typename eT >
inline
grid3D< complex< eT > >::grid3D(const size_t& rows, const size_t& cols, const size_t& lays, const eT& initial)
    : rows(rows)
    , cols(cols)
    , lays(lays)
{
    mem = new complex< eT >[rows * cols * lays];
    complex< eT > init(initial, 0);
    
    std::fill(mem, mem + rows * cols * lays, init);
}

template< typename eT >
inline
grid3D< complex< eT > >::grid3D(const size_t& rcl, const complex< eT >& initial)
    : rows(rcl)
    , cols(rcl)
    , lays(rcl)
{
    mem = new complex< eT >[rows * cols * lays];
    
    std::fill(mem, mem + rows * cols * lays, initial);
}

template< typename eT >
inline
grid3D< complex< eT > >::grid3D(const size_t& rcl, const eT& initial)
    : rows(rcl)
    , cols(rcl)
    , lays(rcl)
{
    mem = new complex< eT >[rows * cols * lays];
    complex< eT > init(initial, 0);
    
    std::fill(mem, mem + rows * cols * lays, init);
}

template< typename eT >
inline
grid3D< complex< eT > >::grid3D(const grid3D< eT >& c)
    : rows(c.rows)
    , cols(c.cols)
    , lays(c.lays)
{
    mem = new complex< eT >[rows * cols * lays];
    
    size_t i;
    for (i = 0; i < rows * cols * lays; ++i)
    {
        mem[i] = complex< eT >(c.mem[i], 0);
    }
}

template< typename eT >
inline
grid3D< complex< eT > >::grid3D(const grid3D< complex< eT > >& c)
    : rows(c.rows)
    , cols(c.cols)
    , lays(c.lays)
{
    mem = new complex< eT >[rows * cols * lays];
    memcpy(mem, c.mem, rows * cols * lays * sizeof(complex< eT >));
}

template< typename eT >
inline
grid3D< complex< eT > >::grid3D(grid3D< complex< eT > >&& c)
    : rows(c.rows)
    , cols(c.cols)
    , lays(c.lays)
{
    complex< eT >* tmp = mem;
    mem                = c.mem;
    c.mem              = tmp;
}

template< typename eT >
inline
grid3D< complex< eT > >::~grid3D()
{
    delete [] mem;
}

template< typename eT >
inline
const grid3D< complex< eT > >& grid3D< complex< eT > >::operator=(const grid3D< eT >& c)
{
    rows = c.rows;
    cols = c.cols;
    lays = c.lays;
    
    delete [] mem;
    mem = new complex< eT >[rows * cols * lays];
    
    size_t i;
    for (i = 0; i < rows * cols * lays; ++i)
    {
        mem[i] = complex< eT >(c.mem[i], 0);
    }
}

template< typename eT >
inline
const grid3D< complex< eT > >& grid3D< complex< eT > >::operator=(const grid3D< complex< eT > >& c)
{
    rows = c.rows;
    cols = c.cols;
    lays = c.lays;
    
    delete [] mem;
    mem = new complex< eT >[rows * cols * lays];
    
    memcpy(mem, c.mem, rows * cols * lays * sizeof(complex< eT >));
}

template< typename eT >
inline
const grid3D< complex< eT > >& grid3D< complex< eT > >::operator=(grid3D< complex< eT > >&& c)
{
    rows = c.rows;
    cols = c.cols;
    lays = c.lays;
    
    complex< eT >* tmp = mem;
    mem                = c.mem;
    c.mem              = tmp;
}

template< typename eT >
inline
complex< eT >& grid3D< complex< eT > >::operator()(const size_t& row, const size_t& col, const size_t& lay)
{
    return mem[rows * cols * lay + cols * col + row];
}

template< typename eT >
inline
const complex< eT >& grid3D< complex< eT > >::operator()(const size_t& row, const size_t& col, const size_t& lay) const
{
    return mem[rows * cols * lay + cols * col + row];
}

template< typename eT >
inline
constexpr size_t grid3D< complex< eT > >::n_rows() const
{
    return rows;
}

template< typename eT >
inline
constexpr size_t grid3D< complex< eT > >::n_cols() const
{
    return cols;
}

template< typename eT >
inline
constexpr size_t grid3D< complex< eT > >::n_lays() const
{
    return lays;
}

template< typename eT >
inline
void grid3D< complex< eT > >::layer_wise_DFT2(const complex< double >& scale)
{
    // declare variables
    size_t i;
    double* data;
    
    // get correct data
    if (is_double< eT >::value == true)
    {
        // If the POD type is double we can just cast the
        // complex array to an double because memory layout
        // is guaranteed by the compiler
        data = reinterpret_cast< double* >(mem);
    }
    else
    {
        // If the POD type is not double we have to create a
        // copy of the memory and cast each element directly
        data = new double[2 * rows * cols * lays];
        
        // Every second element is real or complex. Extract
        // them and store them in data.
        for (i = 0; i < rows * cols * lays; ++i)
        {
            data[i * 2]     = static_cast< double >(mem[i].re);
            data[i * 2 + 1] = static_cast< double >(mem[i].im);
        }
    }
    
    // perform layerwise FFT2
    uzl_fftw_layer_wise_DFT2_grid3D(cols, rows, lays, data);
    
    // skip if scale is default
    if (scale.re != 1 || scale.im != 0)
    {
        // scale data
        if (is_double< eT >::value == true)
        {
            for (i = 0; i < rows * cols * lays; ++i)
            {
                mem[i] *= scale;
            }
        }
        else
        {
            for (i = 0; i < rows * cols * lays; ++i)
            {
                mem[i] = complex< eT >(data[i * 2], data[i * 2 + 1]) * scale;
            }
        }
    }
    
    // free allocated memory
    if (is_double< eT >::value == false)
    {
        delete [] data;
    }
}

template< typename eT >
inline
void grid3D< complex< eT > >::layer_wise_IDFT2(const complex< double >& scale)
{
    // declare variables
    size_t i;
    double* data;
    
    // get correct data
    if (is_double< eT >::value == true)
    {
        // If the POD type is double we can just cast the
        // complex array to an double because memory layout
        // is guaranteed by the compiler
        data = reinterpret_cast< double* >(mem);
    }
    else
    {
        // If the POD type is not double we have to create a
        // copy of the memory and cast each element directly
        data = new double[2 * rows * cols * lays];
        
        // Every second element is real or complex. Extract
        // them and store them in data.
        for (i = 0; i < rows * cols * lays; ++i)
        {
            data[i * 2]     = static_cast< double >(mem[i].re);
            data[i * 2 + 1] = static_cast< double >(mem[i].im);
        }
    }
    
    // perform layerwise FFT2
    uzl_fftw_layer_wise_IDFT2_grid3D(cols, rows, lays, data);
    
    // skip if scale is default
    if (scale.re != 1 || scale.im != 0)
    {
        // scale data
        if (is_double< eT >::value == true)
        {
            for (i = 0; i < rows * cols * lays; ++i)
            {
                mem[i] *= scale;
            }
        }
        else
        {
            for (i = 0; i < rows * cols * lays; ++i)
            {
                mem[i] = complex< eT >(data[i * 2], data[i * 2 + 1]) * scale;
            }
        }
    }
    
    // free allocated memory
    if (is_double< eT >::value == false)
    {
        delete [] data;
    }
}

template< typename eT >
inline
complex< eT >* grid3D< complex< eT > >::memptr()
{
    return mem;
}

template< typename eT >
inline
const complex< eT >* grid3D< complex< eT > >::memptr() const
{
    return mem;
}

template< typename S >
std::ostream& operator<<(std::ostream& o, const grid3D< complex< S > >& c)
{
    std::ios::fmtflags f( std::cout.flags() );
    o << std::endl;
    
    int width   = 20;
    auto format = std::fixed;
    
    if (is_float< S >::value == false && is_double< S >::value == false && is_ldouble< S >::value == false)
    {
        width = 10;
    }
    
    // check values
    size_t x, y, z;
    for (z = 0; z < c.n_lays(); ++z)
    {
        for (x = 0; x < c.n_rows(); ++x)
        {
            for (y = 0; y < c.n_cols(); ++y)
            {
                complex< S > val = c(x, y, z);
                if (std::abs(val.re) >= 10 || std::abs(val.im) >= 10)
                {
                    width   = 22;
                    format  = std::fixed;
                    
                    if (is_float< S >::value == false && is_double< S >::value == false && is_ldouble< S >::value == false)
                    {
                        width = 12;
                    }
                }
                
                if (std::abs(val.re) >= 100 || std::abs(val.im) >= 100)
                {
                    width   = 24;
                    format  = std::fixed;
                    
                    if (is_float< S >::value == false && is_double< S >::value == false && is_ldouble< S >::value == false)
                    {
                        width = 14;
                    }
                }
                
                if (std::abs(val.re) >= 1000 || std::abs(val.im) >= 1000)
                {
                    width   = 28;
                    format  = std::scientific;
                    
                    if (is_float< S >::value == false && is_double< S >::value == false && is_ldouble< S >::value == false)
                    {
                        width = 18;
                    }
                }
            }
        }
    }
    
    // setting decimal precesion
    for (z = 0; z < c.n_lays(); ++z)
    {
        // print layer number
        o << "layer[" << z << "]" << std::endl;
        
        // print numbers of layer
        for (x = 0; x < c.n_rows(); ++x)
        {
            for (y = 0; y < c.n_cols(); ++y)
            {
                // get entry
                complex< S > val = c(x, y, z);
                
                // create string
                std::ostringstream out;
                
                // add real value to string
                out << format << std::setprecision(4) << val.re;
                out << (val.im < 0 ? " - " : " + ") << (val.im == 0 ?  0 : std::abs(val.im)) << "i";
                
                // get string from steram
                std::string str = out.str();
                
                // set filling character
                o << std::setfill(' ') << std::right << std::setw(width) << str;
                
            }
            o << std::endl;
        }
        o << std::endl;
    }
    
    std::cout.flags( f );
    return o;
}

#endif
