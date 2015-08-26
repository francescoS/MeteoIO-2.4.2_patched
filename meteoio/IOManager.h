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
#ifndef __IOMANAGER_H__
#define __IOMANAGER_H__

#include <meteoio/DataGenerator.h>
#include <meteoio/Meteo2DInterpolator.h>
#include <meteoio/BufferedIOHandler.h>
#include <meteoio/MeteoProcessor.h>
#include <meteoio/MeteoData.h>
#include <meteoio/Coords.h>

namespace mio {

class IOManager {

	public:
		enum ProcessingLevel {
			raw           = 1,
			filtered      = 1 << 1,
			resampled     = 1 << 2,
			generated     = 1 << 3,
			num_of_levels = 1 << 4
		};

		IOManager(const std::string& filename_in);
		IOManager(const Config& i_cfg);

		//Legacy support to support functionality of the IOInterface superclass:
		void read2DGrid(Grid2DObject& grid_out, const std::string& parameter="");
		void read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);
		void readDEM(DEMObject& dem_out);
		void readAssimilationData(const Date& date_in, Grid2DObject& da_out);
		void readLanduse(Grid2DObject& landuse_out);
		void readPOI(std::vector<Coords>& pts);
		void write2DGrid(const Grid2DObject& grid_in, const std::string& options="");
		void write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date);
		//end legacy support

		size_t getStationData(const Date& date, STATIONS_SET& vecStation);

		/**
		* @brief Fill vecMeteo with a time series of objects
		* corresponding to the interval indicated by dateStart and dateEnd.
		* Depending on the ProcessingLevel for the instance of the IOManager
		* the data returned will be either raw (read directly from the IOHandler)
		* or processed (read from an BufferedIOHandler and filtered through the
		* MeteoProcessor
		*
		* vecMeteo will be empty if no datasets were retrieved in the interval defined
		* by dateStart and dateEnd
		*
		* Example Usage:
		* @code
		* vector< vector<MeteoData> > vecMeteo;      //empty vector
		* Date d1(2008,06,21,11,00);       //21.6.2008 11:00
		* Date d2(2008,07,21,11,00);       //21.7.2008 11:00
		* IOManager iom(Config("io.ini"));
		* unsigned int nstations = iom.getMeteoData(d1, d2, vecMeteo);
		* @endcode
		* @param dateStart   A Date object representing the beginning of an interval (inclusive)
		* @param dateEnd     A Date object representing the end of an interval (inclusive)
		* @param vecVecMeteo A vector of vector<MeteoData> objects to be filled with data
		* @return            Number of stations for which data has been found in the interval
		*/
		size_t getMeteoData(const Date& dateStart, const Date& dateEnd,
		                    std::vector< METEO_SET >& vecVecMeteo);

		/**
		 * @brief Fill vector<MeteoData> object with multiple instances of MeteoData
		 * corresponding to the instant indicated by a Date object. Each MeteoData
		 * instance within the vector represents the data for one station at the given
		 * instant. Depending on the ProcessingLevel configured data will be either
		 * raw (read directly from the IOHandler)
		 *
		 * NOTE:
		 * - vecMeteo will be empty if there is no data found for any station
		 *
		 * Example Usage:
		 * @code
		 * vector<MeteoData> vecMeteo;      //empty vector
		 * IOManager iomanager(Config("io.ini"));
		 * iomanager.getMeteoData(Date(2008,06,21,11,00), vecMeteo); //21.6.2008 11:00
		 * @endcode
		 * @param i_date      A Date object representing the date/time for the sought MeteoData objects
		 * @param vecMeteo    A vector of MeteoData objects to be filled with data
		 * @return            Number of stations for which data has been found in the interval
		 */
		size_t getMeteoData(const Date& i_date, METEO_SET& vecMeteo);

		/**
		 * @brief Push a vector of time series of MeteoData objects into the IOManager. This overwrites
		 *        any internal buffers that are used and subsequent calls to getMeteoData or interpolate
		 *        will be performed upon this data. This method is a way to bypass the internal reading
		 *        of MeteoData from a certain source and is useful in case the user is only interested
		 *        in data processing and interpolation performed by the IOManager object.
		 * @param level Level of processing that has already been performed on the data (raw XOR filtered)
		 * @param date_start Representing the beginning of the data
		 * @param date_end Representing the end of the data
		 * @param vecMeteo The actual data being pushed into the IOManager object
		 */
		void push_meteo_data(const ProcessingLevel& level, const Date& date_start, const Date& date_end,
		                     const std::vector< METEO_SET >& vecMeteo);

		/**
		 * @brief Fill Grid2DObject with spatial data.
		 * Depending on which meteo plugin is in use, this might be spatially interpolated
		 * point measurements or grids as provided by the data source itself.
		 * Depending on the ProcessingLevel configured data will be either
		 * raw (read directly from the IOHandler)
		 *
		 * NOTE:
		 * - grid will be empty if there is no data found
		 *
		 * Example Usage:
		 * @code
		 * Grid2DObject grid;      //empty grid
		 * IOManager iomanager(Config("io.ini"));
		 * iomanager.getMeteoData(Date(2008,06,21,11,00), MeteoData::TA, grid); //21.6.2008 11:00
		 * @endcode
		 * @param date A Date object representing the date/time for the sought MeteoData objects
		 * @param dem Digital Elevation Model data
		 * @param meteoparam which meteo parameter to return
		 * @param result grid returned filled with the requested data
		 * @return true if the grid got filled
		 */

		bool getMeteoData(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
		                 Grid2DObject& result);

		bool getMeteoData(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
		                 Grid2DObject& result, std::string& info_string);

		void interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
				 const std::vector<Coords>& in_coords, std::vector<double>& result);

		void interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
				 const std::vector<Coords>& in_coords, std::vector<double>& result, std::string& info_string);

		/**
		 * @brief Set the desired ProcessingLevel of the IOManager instance
		 *        The processing level affects the way meteo data is read and processed
		 *        Three values are possible:
		 *        - IOManager::raw data shall be read directly from the buffer
		 *        - IOManager::filtered data shall be filtered before returned to the user
		 *        - IOManager::resampled data shall be resampled before returned to the user
		 *          this only affects the function getMeteoData(const Date&, METEO_DATASET&);
		 *
		 *        The three values can be combined: e.g. IOManager::filtered | IOManager:resampled
		 * @param i_level The ProcessingLevel values that shall be used to process data
		 */
		void setProcessingLevel(const unsigned int& i_level);

		/**
		 * @brief Set buffer window properties requirements as known to the application itself.
		 * This will compare these requirements with the ones expressed by the end user and keep the max between them.
		 * The method can be called several times, it will NOT reset the calculated buffer's requirements but keep
		 * on merging with new submissions. Any parameter given as IOUtils::nodata will be ignored.
		 * @param buffer_size buffer size in days
		 * @param buff_before buffer centering in days
		 */
		void setMinBufferRequirements(const double& buffer_size, const double& buff_before);

		/**
		 * @brief Returns the average sampling rate in the data.
		 * This computes the average sampling rate of the data that is contained in the buffer. This is a quick
		 * estimate, centered on how often a station measures "something" (ie, how many timestamps do we have
		 * for this station in the buffer). if the station measures TA at h+0 and h+30 and
		 * RH at h+15 and h+45, it would return 4 measurements per hour. If the station measures TA and RH at h+0 and h+30,
		 * it would return 2 measurements per hour.
		 * @return average sampling rate in Hz, nodata if the buffer is empty
		 */
		double getAvgSamplingRate() const;

		void writeMeteoData(const std::vector< METEO_SET >& vecMeteo, const std::string& name="");

		/**
		 * @brief Returns a copy of the internal Config object.
		 * This is convenient to clone an iomanager
		 * @return new Config object as a copy of the internal Config
		 */
		const Config getConfig() const;

		const std::string toString() const;

		/**
		 * @brief Add a METEO_SET for a specific instance to the point cache. This is a way to manipulate
		 * MeteoData variables and be sure that the manipulated values are later used for requests
		 * regarding that specific date (e.g. 2D interpolations)
		 *
		 * @param i_date Representing a point in time
		 * @param vecMeteo A vector of MeteoData objects to be copied into the point cache
		 */
		void add_to_cache(const Date& i_date, const METEO_SET& vecMeteo);

		/**
		 * @brief Clear the point cache. All resampled values are dismissed, will need to be recalculated.
		 */
		void clear_cache();

	private:
		void initIOManager();
		void initVirtualStations();
		void fill_filtered_cache();
		bool read_filtered_cache(const Date& start_date, const Date& end_date,
		                         std::vector< METEO_SET >& vec_meteo);
		size_t getTrueMeteoData(const Date& i_date, METEO_SET& vecMeteo);
		size_t getVirtualMeteoData(const Date& i_date, METEO_SET& vecMeteo);

		const Config cfg; ///< we keep this Config object as full copy, so the original one can get out of scope/be destroyed
		IOHandler rawio;
		BufferedIOHandler bufferedio;
		MeteoProcessor meteoprocessor;
		Meteo2DInterpolator interpolator;
		DataGenerator dataGenerator;

		std::vector<size_t> v_params; ///< Parameters for virtual stations
		std::vector<Coords> v_coords; ///< Coordinates for virtual stations
		std::vector<StationData> v_stations; ///< metadata for virtual stations

		ProcessingProperties proc_properties; ///< buffer constraints in order to be able to compute the requested values
		std::map<Date, METEO_SET > virtual_point_cache;  ///< stores already resampled virtual data points
		std::map<Date, METEO_SET > point_cache;  ///< stores already resampled data points
		std::vector< METEO_SET > filtered_cache; ///< stores already filtered data intervals
		Date fcache_start, fcache_end; ///< store the beginning and the end date of the filtered_cache
		unsigned int processing_level;
		bool virtual_stations; ///< compute the meteo values at virtual stations
		bool skip_virtual_stations; ///< skip virtual stations in subsequent calls to prevent recursive calls...
		bool interpol_use_full_dem; ///< use full dem for point-wise spatial interpolations
};
} //end namespace
#endif
