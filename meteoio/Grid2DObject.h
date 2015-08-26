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
#ifndef __GRID2DOBJECT_H__
#define __GRID2DOBJECT_H__

#include <meteoio/Coords.h>
#include <meteoio/Array2D.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/IOUtils.h>

#include <iostream>

namespace mio {

/**
 * @class Grid2DObject
 * @brief A class to represent 2D Grids. Typical application as DEM or Landuse Model.
 *
 * @ingroup data_str
 * @author Thomas Egger
 * @date   2008-12-20
 */

#ifdef _POPC_
#include <paroc_base.h>
class Grid2DObjetctDummy {}; //HACK for POPC

class Grid2DObject : POPBase {
	public:
		/**
		* @brief Serialize method for POPC. Used to marshall data and send it from an object to another.
		*
		* @param buf pointer to the communication buffer
		* @param pack indicates if the data is sent or received
		*/
		virtual void Serialize(POPBuffer &buf, bool pack);
#else
class Grid2DObject {
#endif
	public:
		///structure to contain the grid coordinates of a point in a 2D grid
		typedef struct GRID_POINT_2D {
			size_t ix; ///<grid index along X
			size_t iy; ///<grid index along Y
		} grid_point_2d;

		Grid2DObject& operator=(const Grid2DObject&); ///<Assignement operator
		double& operator ()(const size_t& ix, const size_t& iy);
		double operator ()(const size_t& ix, const size_t& iy) const;
		double& operator ()(const size_t& i);
		double operator ()(const size_t& i) const;

		const std::string toString() const;
		friend std::iostream& operator<<(std::iostream& os, const Grid2DObject& grid);
		friend std::iostream& operator>>(std::iostream& is, Grid2DObject& grid);

		/**
		* @brief Default constructor.
		* Initializes all variables to 0, except lat/long which are initialized to IOUtils::nodata
		*/
		Grid2DObject();
		Grid2DObject(const size_t& ncols, const size_t& nrows,
		             const double& cellsize, const Coords& i_llcorner);

		Grid2DObject(const size_t& ncols, const size_t& nrows,
		             const double& cellsize, const Coords& i_llcorner, const double& init);

		Grid2DObject(const Grid2DObject& i_grid, const double& init);

		Grid2DObject(const size_t& ncols, const size_t& nrows,
		             const double& cellsize, const Coords& i_llcorner, const Array2D<double>& grid2D_in);

		virtual ~Grid2DObject() {};

		/**
		* @brief constructs an object as a subset of another grid object
		* @param i_grid2Dobj initial grid object
		* @param i_nx starting column of the subset
		* @param i_ny starting row of the subset
		* @param i_ncols number of columns of the subset
		* @param i_nrows number of rows of the subset
		*/
		Grid2DObject(const Grid2DObject& i_grid2Dobj,
		             const size_t& i_nx, const size_t& i_ny, //Point in the plane
		             const size_t& i_ncols, const size_t& i_nrows); //dimensions of the sub-plane

		/**
		* @brief Compute the positional parameters that are not already known
		* This means that the Coords::point object that is given either contains geographic coordinates or
		* grid indices. This method will calculate the missing ones (so that (i,j) match with (lat,lon)
		* and (east,north)). If the given point had a "NULL" projection, it will be set to the grid's.
		* @param point coordinate to convert
		* @return false if the given point was invalid or outside the grid (sets (i,j) to closest values within the grid)
		*/
		bool gridify(Coords& point) const;

		/**
		* @brief Compute the positional parameters that are not already known
		* This means that the Coords::point object that is given either contains geographic coordinates or
		* grid indices. This method will calculate the missing ones (so that (i,j) match with (lat,lon)
		* and (east,north)). Any point that is either invalid or outside the grid is removed from the vector.
		* If the given point had a "NULL" projection, it will be set to the grid's.
		* @param vec_points vector containing the coordinates to convert
		* @return false if invalid or external points had to be removed
		*/
		bool gridify(std::vector<Coords>& vec_points) const;

		/**
		* @brief Set all variables in one go.
		* @param ncols number of colums in the grid2D
		* @param nrows number of rows in the grid2D
		* @param cellsize value for cellsize in grid2D
		* @param i_llcorner lower left corner point
		*/
		void set(const size_t& ncols, const size_t& nrows,
		         const double& cellsize, const Coords& i_llcorner);
		/**
		* @brief Set all variables in one go. Notably the member grid2D of type
		* Array2D\<double\> will be destroyed and recreated to size ncols x nrows.
		* @param ncols number of colums in the grid2D
		* @param nrows number of rows in the grid2D
		* @param cellsize value for cellsize in grid2D
		* @param i_llcorner lower left corner point
		* @param grid2D_in grid to be copied by value
		*/
		void set(const size_t& ncols, const size_t& nrows,
		         const double& cellsize, const Coords& i_llcorner, const Array2D<double>& grid2D_in); //TODO: const CArray would be better...

		void set(const size_t& ncols, const size_t& nrows,
		         const double& cellsize, const Coords& i_llcorner, const double& init);

		void size(size_t& o_ncols, size_t& o_nrows) const;
		size_t getNx() const;
		size_t getNy() const;

		/**
		* @brief deletes the data, but keeps geolocalization
		*/
		void clear();

		/**
		* @brief Check if a grid does not contain any data (but it can contain geolocalization)
		* @return true if the grid is 0x0
		*/
		bool isEmpty() const;

		/**
		* @brief check if the current Grid2DObject has the same geolocalization attributes
		* as another Grid2DObject (as well as same cells). The grid coordinates (xllcorner & yllcorner) are NOT
		* checked as these might be tweaked for convenience (like between input grid and local grid)
		* @param target grid to compare to
		* @return true if same geolocalization
		*/
		bool isSameGeolocalization(const Grid2DObject& target) const;

		/**
		* @brief Partitional algorithm to classify each point of the grid.
		* The classification is given by a list of growing thresholds, the 'clusters' are then a simple
		* range of values. Each cluster comes with an 'id' that replaces the values of the points.
		*
		*
		* @param thresholds (const std::vector<double>&) ordered list of thresholds representing a scale of values. Each level of this scale defines a cluster
		* @param ids (const std::vector<double>&) clusters Ids to be used. clustersId.size()=thresholds.size()+1
		* @return true if clusturization was succesfull
		*/
		bool clusterization(const std::vector<double>& thresholds, const std::vector<double>& ids);

		Array2D<double> grid2D; ///<the grid itself (simple 2D table containing the values for each point)
		Coords llcorner; ///<lower left corner of the grid
		double cellsize; ///<dimension in meters of a cell (considered to be square)
		size_t ncols; ///<number of columns in the grid //HACK should point to grid2D.nx
		size_t nrows; ///<number of rows in the grid

 protected:
		void setValues(const size_t& ncols, const size_t& nrows,
		               const double& cellsize, const Coords& i_llcorner);
		void setValues(const size_t& i_ncols, const size_t& i_nrows,
		               const double& i_cellsize);

		/**
		* @brief Converts WGS84 coordinates into grid coordinates (i,j)
		* @param point coordinate to convert
		* @return false if the given point was outside the grid (sets (i,j) to closest values within the grid)
		*/
		bool WGS84_to_grid(Coords& point) const;

		/**
		* @brief Converts grid coordinates (i,j) into WGS84 coordinates
		* @param point coordinate to convert
		* @return false if the given point was invalid (outside the grid or nodata and if possible sets (i,j) to closest values within the grid)
		*/
		bool grid_to_WGS84(Coords& point) const;
};
} //end namespace

#endif
