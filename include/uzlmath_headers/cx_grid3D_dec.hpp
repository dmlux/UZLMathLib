//
//  cx_grid3D_dec.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 08.06.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_cx_grid3D_dec_hpp
#define UZLMathLib_cx_grid3D_dec_hpp

UZLMATH_BEGIN

/*!
 * @ingroup     grid3D
 * @{
 */

/*!
 * @brief       The grid3D class represents a 3D grid containing data of a
 *              given type.
 * @details     The grid3D is a 3D data structure that uses matrices in
 *              several layers.
 *
 * @tparam      T The type of data the grid is holding
 *
 * @since       0.0.1
 *
 * @date        19.05.2015
 * @author      Denis-Michael Lux <denis.lux@icloud.com>
 *
 */
template< typename T >
struct
grid3D< complex< T >, if_pod_type< T > >
{
    // typedefs
    typedef T pod_type;             //!< POD type of the elements
    
    // ivars
    const size_t rows;              //!< number of rows in each layer
    const size_t cols;              //!< number of cols in each layer
    const size_t lays;              //!< number of layers
    
    const complex< pod_type >* mem; //!< storage of the 3D grid
    
    // methods
    inline                                      grid3D();
    inline                                      grid3D(const size_t& rcl);
    inline                                      grid3D(const size_t& rcl, const complex< pod_type >& initial);
    inline                                      grid3D(const size_t& rcl, const pod_type& initial);
    inline                                      grid3D(const size_t& rows, const size_t& cols, const size_t& lays);
    inline                                      grid3D(const size_t& rows, const size_t& cols, const size_t& lays, const complex< pod_type >& initial);
    inline                                      grid3D(const size_t& rows, const size_t& cols, const size_t& lays, const pod_type& initial);
    inline                                      grid3D(const grid3D< pod_type >& c);
    inline                                      grid3D(const grid3D< complex< pod_type > >& c);
    inline                                      grid3D(grid3D< complex< pod_type > >&& c);
    inline                                     ~grid3D();
    
    inline const grid3D< complex< pod_type > >& operator=(const grid3D< pod_type >& c);
    inline const grid3D< complex< pod_type > >& operator=(const grid3D< complex< pod_type > >& c);
    inline const grid3D< complex< pod_type > >& operator=(grid3D< complex< pod_type > >&& c);
    
    inline       complex< pod_type >&           operator()(const size_t& row, const size_t& col, const size_t& lay);
    inline const complex< pod_type >&           operator()(const size_t& row, const size_t& col, const size_t& lay) const;
    
    inline       void                           layer_wise_DFT2(const complex< double >& scale = complex< pod_type >(1, 0));
    inline       void                           layer_wise_IDFT2(const complex< double >& scale = complex< pod_type >(1, 0));
};

template< typename S > std::ostream& operator<<(std::ostream& o, const grid3D< complex< S > >& c);

/*!
 * @}
 */

UZLMATH_END

#endif /* cx_grid3D_dec.hpp */
