/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef DEMOBJECT_H
#define DEMOBJECT_H

#include <meteoio/Array2D.h>
#include <meteoio/Grid2DObject.h>
#include <meteoio/IOUtils.h>
#include <meteoio/meteolaws/Meteoconst.h> //for math constants

#include <cmath>
#include <limits>

namespace mio {

/**
 * @class DEMObject
 * @brief A class to represent DEMs and automatically compute some properties.
 * This class stores elevation grids and their georeferencing, expressed as the lower-left coordinates, the cellsize (cells are assumed to be square) and a nodata code for potentially empty cells (The nodata parameter is supposed to be IOUtils::nodata).
 * This class also automatically computes local slope, azimuth, curvature, normals and minimal/maximal for normalization.
 * Various algorithms are available to compute these properties (see mio::DEMObject::slope_type) and it is possible to toggle between automatic refresh or not. Several other DEM related values can be computed, such as the horizon, displacements within the DEM, etc
 *
 * @ingroup data_str
 * @author Gaël Rosset - Mathias Bavay
 * @date   2009-07-20
 */

#ifdef _POPC_
class DEMObjectDummy {}; //HACK for POPC

#include <paroc_base.h>
class DEMObject : public Grid2DObject/*, POPBase*/ {
	public:
		void Serialize(POPBuffer &buf, bool pack);
#else
class DEMObject : public Grid2DObject {
#endif
	public:
		Array2D<double> slope;
		Array2D<double> azi;
		Array2D<double> curvature;
		Array2D<double> Nx, Ny, Nz;
		double min_altitude, min_slope, min_curvature;
		double max_altitude, max_slope, max_curvature;

		///Keywords for slope computation algorithm
		typedef enum SLOPE_TYPE {
			DFLT, ///< whatever algorithm that has been defined as default
			FLEM, ///< four nearest neighbors (Fleming and Hoffer, 1979). It seems to be the same as (Zevenbergen and Thorne, 1987)
			HICK, ///< maximum downhill slope method (Dunn and Hickey, 1998)
			HORN, ///< eight neighbor algorithm (Horn, 1981) as used by ArcGIS. It seems to be the same as (Corripio, 2002) but for border cells.
			CORR, ///< surface normal vector using the two triangle method (Corripio, 2003) and eight-neighbor algorithm (Horn, 1981) for border cells
			D8 ///< discretized azimuth directions (angles for N, NE, etc) and slope rounded to nearest integer
		} slope_type;

		///Keywords for automatic update of parameters. They can be combined with "|"
		typedef enum UPDATE_TYPE {
			NO_UPDATE=0, ///< no updates at all
			SLOPE=1, ///< update the slopes
			NORMAL=2, ///< update the normals
			CURVATURE=4 ///< update the curvatures
		} update_type;

		DEMObject(const slope_type& i_algorithm=DFLT);

		DEMObject(const size_t& ncols_in, const size_t& nrows_in,
		          const double& cellsize_in, const Coords& llcorner_in, const slope_type& i_algorithm=DFLT);

		DEMObject(const size_t& ncols_in, const size_t& nrows_in,
		          const double& cellsize_in, const Coords& llcorner_in, const Array2D<double>& altitude_in,
		          const bool& i_update=true, const slope_type& i_algorithm=DFLT);

		DEMObject(const Grid2DObject& dem_in, const bool& i_update=true, const slope_type& i_algorithm=DFLT);

		DEMObject (const DEMObject& i_dem,
		           const size_t& i_nx, const size_t& i_ny, //Point in the plane
		           const size_t& i_ncols, const size_t& i_nrows, //dimensions of the sub-plane
		           const bool& i_update=true, const slope_type& i_algorithm=DFLT);

		void setDefaultAlgorithm(const slope_type& i_algorithm);
		int getDefaultAlgorithm() const;
		void setUpdatePpt(const update_type& in_update_flag);
		int getUpdatePpt() const;

		void update(const std::string& algorithm);
		void update(const slope_type& algorithm=DFLT);
		void updateAllMinMax();
		void printFailures();
		void sanitize();

		Grid2DObject getHillshade(const double& elev=38., const double& azimuth=0.) const;
		double horizontalDistance(const double& xcoord1, const double& ycoord1, const double& xcoord2, const double& ycoord2);
		double horizontalDistance(Coords point1, const Coords& point2);
		double terrainDistance(Coords point1, const Coords& point2);
		void getPointsBetween(Coords point1, Coords point2, std::vector<GRID_POINT_2D>& vec_points);
		void getPointsBetween(const Coords& point, const double& bearing, std::vector<GRID_POINT_2D>& vec_points);
		double getHorizon(const Coords& point, const double& bearing);
		void getHorizon(const Coords& point, const double& increment, std::vector<double>& horizon);

		friend std::iostream& operator<<(std::iostream& os, const DEMObject& dem);
		friend std::iostream& operator>>(std::iostream& is, DEMObject& dem);

	private:
		void CalculateAziSlopeCurve(slope_type algorithm);
		double CalculateAspect(const double& o_Nx, const double& o_Ny, const double& o_Nz, const double& o_slope, const double no_slope=Cst::PI);
		void CalculateHick(double A[4][4], double& o_slope, double& o_Nx, double& o_Ny, double& o_Nz);
		void CalculateFleming(double A[4][4], double& o_slope, double& o_Nx, double& o_Ny, double& o_Nz);
		void CalculateHorn(double A[4][4], double& o_slope, double& o_Nx, double& o_Ny, double& o_Nz);
		void CalculateCorripio(double A[4][4], double& o_slope, double& o_Nx, double& o_Ny, double& o_Nz);
		void (DEMObject::*CalculateSlope)(double A[4][4], double& o_slope, double& o_Nx, double& o_Ny, double& o_Nz);
		double getCurvature(double A[4][4]);

		double steepestGradient(double A[4][4]);
		double lineGradient(const double& A1, const double& A2, const double& A3);
		double fillMissingGradient(const double& delta1, const double& delta2);
		void surfaceGradient(double& dx_sum, double& dy_sum, double A[4][4]);
		double avgHeight(const double& z1, const double &z2, const double& z3);
		void getNeighbours(const size_t i, const size_t j, double A[4][4]);
		double safeGet(const int i, const int j);

		int update_flag;
		slope_type dflt_algorithm;
		size_t slope_failures; ///<contains the number of points that have an elevation but no slope
		size_t curvature_failures; ///<contains the number of points that have an elevation but no curvature
};
} //end namespace

#endif
