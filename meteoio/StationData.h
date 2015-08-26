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
#ifndef __STATIONDATA_H__
#define __STATIONDATA_H__

#include <meteoio/Coords.h>
#include <meteoio/IOUtils.h>

#include <string>
#include <sstream>
#include <iomanip>
#include <vector>

namespace mio {

/**
 * @class StationData
 * @brief A class to represent meteo stations with attributes like longitude, latitude, etc.
 *
 * @ingroup data_str
 * @author Thomas Egger
 * @date   2008-11-29
 */

#ifdef _POPC_
#include <paroc_base.h>
class StationData :POPBase {
	public:
		void Serialize(POPBuffer &buf, bool pack);
#else
class StationData {
#endif
	public:
		//Constructors
		/**
		* @brief The default constructor initializing every double attribute to nodata and strings to  ""
		*/
		StationData(void);

		/**
		* @brief A constructor that takes three to six arguments
		* @param i_position Position of the station
		* @param i_id Station's id (short identification)
		* @param i_name Full name of the station (default "")
		*/
		StationData(const Coords& i_position, const std::string& i_id="", const std::string& i_name="");

		//Specific getter functions
		std::string getStationID() const;
		std::string getStationName() const;
		Coords getPosition() const;

		/**
		* @brief Get local slope angle
		* @return angle of the local slope (in degrees, between 0 and 90 degrees)
		*/
		double getSlopeAngle() const;

		/**
		* @brief Get local slope azimuth
		* @return azimuth of the local slope expressed as a bearing (0 is North, in degrees, clockwise)
		*/
		double getAzimuth() const;

		/**
		* @brief General setter function
		* @param i_position Position of the station
		* @param i_id Station's id (short identification)
		* @param i_name Name of the station (default "")
		*/
		void setStationData(const Coords& i_position, const std::string& i_id="", const std::string& i_name="");

		/**
		* @brief Slope information setter
		* @param in_slope_angle angle of the local slope (in degrees, between 0 and 90 degrees)
		* @param in_azimuth azimuth of the local slope expressed as a bearing (0 is North, in degrees, clockwise)
		*/
		void setSlope(const double& in_slope_angle, const double& in_azimuth);

		const std::string toString() const;
		friend std::iostream& operator<<(std::iostream& os, const StationData& station);
		friend std::iostream& operator>>(std::iostream& is, StationData& station);

		//Comparison operators
		/**
		* @brief Equality %operator
		* check all parameters but the station name
		* @return true or false
		*/
		bool operator==(const StationData&) const;
		bool operator!=(const StationData&) const; ///<Operator that tests for inequality

		/**
		* @brief Simple merge strategy.
		* If some fields of the first argument are empty, they will be filled by the matching field from the
		* second argument.
		* @param sd1 first StationData to merge, highest priority
		* @param sd2 second StationData to merge, lowest priority
		* @return new StationData object
		*/
		static StationData merge(const StationData& sd1, const StationData& sd2);

		/**
		* @brief Simple merge strategy.
		* If some fields of the current object are empty, they will be filled by the macthing field from the
		* provided argument.
		* @param sd2 extra StationData to merge, lowest priority
		*/
		void merge(const StationData& sd2);

	public:
		Coords position;
		std::string stationID; ///<ID of the Station, typically a short string
		std::string stationName; ///<Name of the Station, detailled description

	private:
		double slope; ///<Local slope at the station, in degrees, between 0 and 90 degrees
		double azi; ///<Azimuth at the local slope at the station, in degrees, 0 at north, compass orientation
};

typedef std::vector<StationData> STATIONS_SET;

} //end namespace

#endif
