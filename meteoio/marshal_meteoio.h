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
#ifdef _POPC_

#ifndef MARSHAL_METEOIO_H
#define MARSHAL_METEOIO_H

#include <meteoio/Grid2DObject.h>
#include <meteoio/Grid3DObject.h>
#include <meteoio/DEMObject.h>
#include <meteoio/StationData.h>
#include <meteoio/MeteoData.h>
#include <meteoio/Coords.h>
#include <meteoio/Config.h>

#include <vector>
#include <string>

#include <paroc_memspool.h>
#include <paroc_buffer.h>

namespace mio {

typedef Array2D<double> DOUBLE2D; //HACK for POPC
typedef Array3D<double> DOUBLE3D;
typedef Array2D<int> INT2D;
typedef Array2D<char> CHAR2D;
typedef std::vector<std::string> STR_VECTOR;

void marshal_uint(POPBuffer &buf,unsigned int &data, int maxsize, int flag, POPMemspool *temp);

void marshal_MeteoParameters(POPBuffer &buf, MeteoData::Parameters &data, int maxsize, int flag, POPMemspool *temp);
void marshal_MeteoGridsParameters(POPBuffer &buf, MeteoGrids::Parameters &data, int maxsize, int flag, POPMemspool *temp);

void marshal_slope_type(POPBuffer &buf, DEMObject::slope_type &data, int maxsize, int flag, POPMemspool *temp);

void marshal_geo_distances(POPBuffer &buf, Coords::geo_distances &data, int maxsize, int flag, POPMemspool *temp);

void marshal_DOUBLE2D(POPBuffer &buf, DOUBLE2D &data,int maxsize, int flag, POPMemspool *temp);

void marshal_DOUBLE3D(POPBuffer &buf, DOUBLE3D &data,int maxsize, int flag, POPMemspool *temp);

void marshal_INT2D(POPBuffer &buf, INT2D &data,int maxsize, int flag, POPMemspool *temp);

void marshal_CHAR2D(POPBuffer &buf, CHAR2D &data,int maxsize, int flag, POPMemspool *temp);

void marshal_vec_coords(POPBuffer &buf,std::vector<Coords> &data, int maxsize, int flag, POPMemspool *temp);

void marshal_METEO_SET(POPBuffer &buf, METEO_SET &data, int maxsize, int flag, POPMemspool *temp);

void marshal_vector_METEO_SET(POPBuffer &buf, std::vector<METEO_SET> &data, int maxsize, int flag, POPMemspool *temp);

void marshal_vecstr(POPBuffer &buf, std::vector<std::string> &data, int maxsize, int flag, POPMemspool *temp);

void marshal_map_str_str(POPBuffer &buf, std::map<std::string, std::string> &data_map, int maxsize, int flag, POPMemspool *temp);

void marshal_map_str_dbl(POPBuffer &buf, std::map<std::string, double> &data_map, int maxsize, int flag, POPMemspool *temp);

void marshal_map_str_vecstr(POPBuffer &buf, std::map<std::string, STR_VECTOR> &data_map, int maxsize, int flag, POPMemspool *temp);

void marshal_Coords(POPBuffer &buf, Coords &data, int maxsize, int flag, POPMemspool *temp);

void marshal_STATIONS_SET(POPBuffer &buf, STATIONS_SET &data, int maxsize, int flag, POPMemspool *temp);

void marshal_vector_STATIONS_SET(POPBuffer &buf, std::vector<STATIONS_SET> &data, int maxsize, int flag, POPMemspool *temp);

void marshal_vector_Grid2DObject(POPBuffer &buf, std::vector<Grid2DObject> &data, int maxsize, int flag, POPMemspool *temp);

void marshal_DEMObject(POPBuffer &buf, DEMObject &data, int maxsize, int flag, POPMemspool *temp);

void marshal_Date(POPBuffer &buf, Date &data, int maxsize, int flag, POPMemspool *temp);

void marshal_Config(POPBuffer &buf, Config &data, int maxsize, int flag, POPMemspool *temp);

void marshal_Grid2DObject(POPBuffer &buf, Grid2DObject &data, int maxsize, int flag, POPMemspool *temp);

void marshal_MeteoData(POPBuffer &buf, MeteoData &data, int maxsize, int flag, POPMemspool *temp);

} //end namespace

#endif

#endif
