//
//  grid3D_dec.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 19.05.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_grid3D_dec_hpp
#define UZLMathLib_grid3D_dec_hpp

UZLMATH_BEGIN

/*!
 * @brief       A collection of classes and functions for 3D grids
 * @details     A 3D grid can hold values of the given template type. The
 * @defgroup    grid3D 3D grid
 * @{
 */

/*!
 * @brief       The grid3D class represents a 3D grid containing data of a
 *              given type.
 * @details     The grid3D is a 3D data structure that uses matrices in
 *              several layers.
 *
 * @tparam      eT The type of data the grid is holding
 *
 * @since       0.0.1
 * 
 * @date        19.05.2015
 * @author      Denis-Michael Lux <denis.lux@icloud.com>
 *
 */
template< typename eT >
struct grid3D
{
    const size_t rows;         //!< number of rows in each layer
    const size_t cols;         //!< number of cols in each layer
    const size_t lays;         //!< number of layers
    
    const eT* mem;             //!< storage of the 3D grid
    
    inline                     grid3D();
    inline                     grid3D(const size_t& rcl);
    inline                     grid3D(const size_t& rcl, const eT& initial);
    inline                     grid3D(const size_t& rows, const size_t& cols, const size_t& lays);
    inline                     grid3D(const size_t& rows, const size_t& cols, const size_t& lays, const eT& initial);
    inline                     grid3D(const grid3D< eT >& c);
    inline                     grid3D(grid3D< eT >&& c);
    inline                    ~grid3D();
    
    inline const grid3D< eT >& operator=(const grid3D< eT >& c);
    inline const grid3D< eT >& operator=(grid3D< eT >&& c);
    
    inline       eT&           operator()(const size_t& row, const size_t& col, const size_t& lay);
    inline const eT&           operator()(const size_t& row, const size_t& col, const size_t& lay) const;
};

template< typename S > std::ostream& operator<<(std::ostream& o, const grid3D< S >& c);

/*!
 * @}
 */

UZLMATH_END

#endif
