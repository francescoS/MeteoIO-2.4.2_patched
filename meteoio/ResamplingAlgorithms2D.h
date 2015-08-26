/***********************************************************************************/
/*  Copyright 2011 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef __RESAMPLINGALGORITHMS2D_H__
#define __RESAMPLINGALGORITHMS2D_H__

#include <meteoio/Grid2DObject.h>
#include <iostream>
#include <string>

namespace mio {

/**
 * @class ResamplingAlgorithms2D
 * @brief Spatial resampling algorithms
 *
 * @ingroup stats
 * @author Mathias Bavay
 * @date   2011-06-29
 */
class ResamplingAlgorithms2D {
	public:
		//Available algorithms
		static const Grid2DObject NearestNeighbour(const Grid2DObject &i_grid, const double &factor);
		static const Grid2DObject BilinearResampling(const Grid2DObject &i_grid, const double &factor);
		static const Grid2DObject cubicBSplineResampling(const Grid2DObject &i_grid, const double &factor);

		static const Array2D<double> NearestNeighbour(const Array2D<double> &i_grid, const double &factor_x, const double &factor_y);
		static const Array2D<double> BilinearResampling(const Array2D<double> &i_grid, const double &factor_x, const double &factor_y);
		static const Array2D<double> cubicBSplineResampling(const Array2D<double> &i_grid, const double &factor_x, const double &factor_y);

	private:
		static void cubicBSpline(Array2D<double> &o_grid, const Array2D<double> &i_grid);
		static void Bilinear(Array2D<double> &o_grid, const Array2D<double> &i_grid);
		static void NearestNeighbour(Array2D<double> &o_grid, const Array2D<double> &i_grid);

		static double bilinear_pixel(const Array2D<double> &i_grid, const size_t &org_ii, const size_t &org_jj, const size_t &org_ncols, const size_t &org_nrows, const double &x, const double &y);
		static double BSpline_weight(const double &x);
};
} //end namespace

#endif
