//
//  fn_obj2file_grid3D.hpp
//  UZLMathLib
//
//  Created by Denis-Michael Lux on 16.06.15.
//
//  This software may be modified and distributed under the terms
//  of the BSD license. See the LICENSE file for details.
//

#ifndef UZLMathLib_fn_obj2file_grid3D_hpp
#define UZLMathLib_fn_obj2file_grid3D_hpp

UZLMATH_BEGIN

/*!
 * @brief           Stores the contents of a complex grid3D object on disk.
 * @details         Since grid objects can be very big this method uses
 *                  the old fashiond C way to store the cube. The resulting
 *                  file contains each element in two rows seperated by a 
 *                  '\n' character. The fastest index is the column index,
 *                  followed by the row index. The slowest index is the layer
 *                  index.
 *
 * @param[in]       grid The complex grid that is supposed to be stored on disk
 * @param[in]       fileName The file name including the path were the 
 *                  file should be stored.
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            16.06.2015
 *
 * @since           0.1.1
 *
 * @ingroup         grid3D
 */
template< typename T >
inline
typename uzl_void_real_num_only< T >::result obj2file(const grid3D< complex< T > >& grid, const std::string& fileName)
{
    // create file and open it with write flag
    FILE* fp = fopen(fileName.c_str(), "w");
    
    // iterate over grid starting with layers
    unsigned int i, j, k;
    for (i = 0; i < grid.n_lays(); ++i)
    {
        
        // from this point iterate over the layer the same way as an matrix.
        // first iterate over the rows of the layer (matrix) and the over the
        // columns.
        for (j = 0; j < grid.n_rows(); ++j)
        {
            for (k = 0; k < grid.n_cols(); ++k)
            {
                fprintf(fp, "%.16f\n%.16f\n", grid(j, k, i).re, grid(j, k, i).im);
            }
        }
    }
    
    // close file and everything is done
    fclose(fp);
}

/*!
 * @brief           Stores the contents of a grid3D object on disk.
 * @details         Since grid objects can be very big this method uses
 *                  the old fashiond C way to store the cube. The resulting
 *                  file contains each element in two rows seperated by a
 *                  '\n' character. The fastest index is the column index,
 *                  followed by the row index. The slowest index is the layer
 *                  index.
 *
 * @param[in]       grid The grid that is supposed to be stored on disk
 * @param[in]       fileName The file name including the path were the
 *                  file should be stored.
 *
 * @author          Denis-Michael Lux <denis.lux@icloud.com>
 * @date            16.06.2015
 *
 * @since           0.1.1
 *
 * @ingroup         grid3D
 */
template< typename T >
inline
typename uzl_void_real_num_only< T >::result obj2file(const grid3D< T >& grid, const std::string& fileName)
{
    // create file and open it with write flag
    FILE* fp = fopen(fileName.c_str(), "w");
    
    // iterate over grid starting with layers
    unsigned int i, j, k;
    for (i = 0; i < grid.n_lays(); ++i)
    {
        
        // from this point iterate over the layer the same way as an matrix.
        // first iterate over the rows of the layer (matrix) and the over the
        // columns.
        for (j = 0; j < grid.n_rows(); ++j)
        {
            for (k = 0; k < grid.n_cols(); ++k)
            {
                fprintf(fp, "%.16f\n", grid(j, k, i));
            }
        }
    }
    
    // close file and everything is done
    fclose(fp);
}

UZLMATH_END

#endif
