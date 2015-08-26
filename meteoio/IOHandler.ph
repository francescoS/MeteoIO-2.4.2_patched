/***********************************************************************************/
/*  Copyright 2009 HES-SO Fribourg                                                 */
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
#ifndef __IOHANDLER_H__
#define __IOHANDLER_H__

#include <meteoio/IOInterface.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/IOPlugin.h>
#include <meteoio/marshal_meteoio.h>

#include <map>
#include <string>

using namespace mio; //HACK for POPC

namespace mio {

typedef std::map<std::string, IOPlugin>::iterator PLUGIN_ITERATOR; //HACK for POPC
//This one line above does absolutely nothing, but if removed, popc does not compile the file....

/**
* @file IOHandler.ph
* The is the parclass implementing the interface as defined by the IOInterface class.
* This class is responsible for loading the necessary plugins and getting the data through them.
*/
parclass IOHandler {
// Note : No heritage here for POPC++ : a parclass cannot herit from a class
		classuid(1003);
	public:
		//IOHandler(const IOHandler&) @{od.url("localhost");};
		IOHandler(const Config&) @{od.url("localhost");}; //@{ power=100 ?: 50; };
		~IOHandler();

		//methods defined in the IOInterface class
		virtual void read2DGrid([out]Grid2DObject& out_grid, const std::string& parameter="");
		virtual void read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);
		virtual void readDEM([out]DEMObject& dem_out);
		virtual void readLanduse([out]Grid2DObject& landuse_out);
		virtual void readStationData([in]const Date& date,
		              [proc=marshal_STATION_TIMESERIE] STATION_TIMESERIE& vecStation);
		virtual void writeMeteoData([in,proc=marshal_vector_METEO_TIMESERIE] std::vector<METEO_TIMESERIE>& vecMeteo,
		              [in]const std::string& name="");
		virtual void readMeteoData([in]const Date& dateStart, [in]const Date& dateEnd,
		              [proc=marshal_vector_METEO_TIMESERIE] std::vector<METEO_TIMESERIE>& vecMeteo,
		              const unsigned& stationindex=IOUtils::npos);
		void readMeteoData([in]const Date& date, [proc=marshal_METEO_TIMESERIE] METEO_TIMESERIE& vecMeteo);
		virtual void readAssimilationData([in] const Date&,[out] Grid2DObject& da_out);
		virtual void readSpecialPoints([out,proc=marshal_vec_coords]std::vector<Coords>& pts);
		virtual void write2DGrid([in]const Grid2DObject& grid_in, [in]const std::string& name);
		virtual void write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date);

		//friend std::ostream& operator<<(std::ostream& os, const IOHandler& data); //not "friends" in a parclass!
		std::string toString()/* const*/; //HACK for POPC

	private:
		void loadDynamicPlugins();
		void loadPlugin(const std::string& libname, const std::string& classname,
		                DynamicLibrary*& dynLibrary, IOInterface*& io);
		void deletePlugin(DynamicLibrary*& dynLibrary, IOInterface*& io);
		void registerPlugins();
		IOInterface *getPlugin(const std::string& cfgkey, const std::string& cfgsection="GENERAL");

		const Config& cfg;
		std::map<std::string, IOPlugin> mapPlugins;

		bool enable_copying;
		std::vector<std::string> copy_parameter, copy_name;
};

} //namespace

#endif
