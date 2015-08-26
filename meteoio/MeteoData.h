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
#ifndef __METEODATA_H__
#define __METEODATA_H__

#include <meteoio/Date.h>
#include <meteoio/StationData.h>
#include <meteoio/IOUtils.h>

#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>

namespace mio {

class MeteoData; //forward declaration
typedef std::vector<MeteoData> METEO_SET;

/**
 * @class MeteoGrids
 * @brief A class to represent the meteorological parameters that could be contained in a grid.
 * This should be very close to MeteoData with a few additions (like the wind u,v,w)
 * @ingroup data_str
 * @author Mathias Bavay
 * @date   2011-12-22
 */

class MeteoGrids {
	public:
		/// \anchor meteogrids this enum provides names for possible meteogrids (from an ARPS file, etc)
		enum Parameters {firstparam=0,
		                 TA=firstparam, ///< Air temperature
		                 RH, ///< Relative humidity
		                 VW, ///< Wind velocity
		                 DW, ///< Wind direction
		                 VW_MAX, ///< Maximum wind velocity
		                 ISWR, ///< Incoming short wave radiation
		                 RSWR, ///< Reflected short wave radiation
		                 ILWR, ///< Incoming long wave radiation
		                 //TAU_CLD,  ///< Cloud transmissivity, or ISWR/ISWR_clear_sky, aka solar index
		                 HS, ///< Height of snow
		                 HNW, ///< Water equivalent of precipitations, either solid or liquid
		                 TSG, ///< Temperature ground surface
		                 TSS, ///< Temperature snow surface
		                 P, ///< Air pressure
		                 U, ///< East component of wind
		                 V, ///< North component of wind
		                 W, ///< Vertical component of wind
		                 SWE, ///< Snow Water Equivalent
		                 ROT, ///< Total generated runoff
		                 ALB, ///< Albedo
		                 DEM, ///< Digital Elevation Model
		                 SLOPE, ///< DEM slope angle
		                 AZI, ///< DEM slope azimuth
		                 lastparam=AZI};

		static const size_t nrOfParameters; ///<holds the number of meteo parameters stored in MeteoData
		static const std::string& getParameterName(const size_t& parindex);

	private:
		//static methods
		static std::vector<std::string> paramname;
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map meteoparamname
};

/**
 * @class MeteoData
 * @brief A class to represent a singular measurement received from one station at a certain time (represented by the Date object)
 *
 * @ingroup data_str
 * @author Thomas Egger
 * @date   2008-12-05
 */

#ifdef _POPC_
#include <paroc_base.h>
class MeteoData : POPBase {
	public:
		void Serialize(POPBuffer &buf, bool pack);
#else
class MeteoData {
#endif
	public:
		/// \anchor meteoparam this enum provides indexed access to meteorological fields
		enum Parameters {firstparam=0,
		                 P=firstparam, ///< Air pressure
		                 TA, ///< Air temperature
		                 RH, ///< Relative humidity
		                 TSG, ///< Temperature of the ground surface
		                 TSS, ///< Temperature of the snow surface
		                 HS, ///< Height of snow
		                 VW, ///< Wind velocity
		                 DW, ///< Wind direction
		                 VW_MAX, ///< Maximum wind velocity
		                 RSWR, ///< Reflected short wave radiation
		                 ISWR, ///< Incoming short wave radiation
		                 ILWR, ///< Incoming long wave radiation (downwelling)
		                 //TAU_CLD,  ///< Cloud transmissivity, or ISWR/ISWR_clear_sky, aka solar index
		                 HNW, ///< Water equivalent (water depth) of precipitations, either solid or liquid
		                 lastparam=HNW};

		static const std::string& getParameterName(const size_t& parindex);

		/**
		 * @brief The default constructor initializing every double attribute to nodata and the Date to julian==0.0
		 */
		MeteoData(void);

		/**
		* @brief A constructor that sets the measurment time
		* @param in_date A Date object representing the time of the measurement
		*/
		MeteoData(const Date& in_date);

		/**
		* @brief A constructor that sets the measurment time and meta data
		* @param date_in A Date object representing the time of the measurement
		* @param meta_in A StationData object containing the meta data
		*/
		MeteoData(const Date& date_in, const StationData& meta_in);

		/**
		* @brief A setter function for the measurement date
		* @param in_date A Date object representing the time of the measurement
		*/
		void setDate(const Date& in_date);

		/**
		* @brief Add another variable to the MeteoData object,
		*        a double value will be added and the nrOfParameters increased
		* @param i_paramname A parameter name, e.g. "VSWR"
		* @return A size_t denoting the index of the the parameter added
		*/
		size_t addParameter(const std::string& i_paramname);

		/**
		* @brief Check whether a certain parameter is a part of this MeteoData instance
		* @param parname A string parameter, representing a meteo parameter, e.g. "VSWR"
		* @return A boolean indicating whether the parameter is a part of the object
		*/
		bool param_exists(const std::string& parname) const;

		/**
		 * @brief Resets all the meteo parameters to IOUtils::nodata
		 *        NOTE: member vars date and resampled are not affected
		 */
		void reset();

		bool isResampled() const;
		void setResampled(const bool&);

		void standardizeNodata(const double& plugin_nodata);

		double& operator()(const size_t& parindex);
		const double& operator()(const size_t& parindex) const;
		double& operator()(const std::string& parname);
		const double& operator()(const std::string& parname) const;

		const std::string& getNameForParameter(const size_t& parindex) const;
		size_t getParameterIndex(const std::string& parname) const;
		size_t getNrOfParameters() const;

		/**
		 * @brief Simple merge strategy for vectors containing meteodata for a given timestamp.
		 * If some fields of the MeteoData objects given in the first vector are nodata, they will be
		 * filled by the matching field from the MeteoData objects given in the second vector (if the
		 * same location exist). Stations only occuring in the second vector will be appended to the
		 * first vector.
		 * @note two stations are considered to be identical if they fit within a 5m 3D box
		 * @note the vectors are supposed to contain data at a given time stamp. If both vectors don't match a
		 * common time stamp, nothing is done
		 * @param vec1 reference vector, highest priority
		 * @param vec2 extra vector to merge, lowest priority
		 * @param simple_merge if set to true, assume all stations are unique (ie. simply append vec2 to vec1)
		 */
		static void merge(std::vector<MeteoData>& vec1, const std::vector<MeteoData>& vec2, const bool& simple_merge=false);

		/**
		 * @brief Simple merge strategy for vectors containing meteodata for a given timestamp.
		 * If some fields of the MeteoData objects given in the first vector are nodata, they will be
		 * filled by the matching field from the MeteoData object given in the second argument (if the
		 * same location exist). If meteo2 does not describe a station already in vec, it will simply be appended.
		 * @note two stations are considered to be identical if they fit within a 5m 3D box
		 * @note the datasets are supposed to contain data at a given time stamp. If vec1 and meteo2 don't match a
		 * common time stamp, nothing is done
		 * @param vec reference vector, highest priority
		 * @param meteo2 extra MeteoData object to merge, lowest priority
		 * @param simple_merge if set to true, assume all stations are unique (ie.simply append meteo2 to vec)
		 */
		static void merge(std::vector<MeteoData>& vec, const MeteoData& meteo2, const bool& simple_merge=false);

		/**
		 * @brief Simple merge strategy.
		 * If some fields of the object given as first argument are nodata, they will be filled by the matching field from the
		 * provided argument.
		 * @note no check on the location is performed, ie. it can merge data from stations kilometers away...
		 * @param meteo1 reference MeteoData, highest priority
		 * @param meteo2 extra MeteoData to merge, lowest priority
		 */
		static MeteoData merge(const MeteoData& meteo1, const MeteoData& meteo2);

		/**
		 * @brief Simple merge strategy.
		 * If some fields of the current object are nodata, they will be filled by the matching field from the
		 * provided argument.
		 * @note no check on the location is performed, ie. it can merge data from stations kilometers away...
		 * @param meteo2 extra MeteoData to merge, lowest priority
		 */
		void merge(const MeteoData& meteo2);

		const std::string toString() const;
		friend std::iostream& operator<<(std::iostream& os, const MeteoData& data);
		friend std::iostream& operator>>(std::iostream& is, MeteoData& data);

		//Comparison operators
		bool operator==(const MeteoData&) const; ///<Operator that tests for equality
		bool operator!=(const MeteoData&) const; ///<Operator that tests for inequality

		//direct access allowed
		Date date; ///<Timestamp of the measurement
		StationData meta; ///<The meta data of the measurement

		static const size_t nrOfParameters; ///<holds the number of meteo parameters stored in MeteoData

	private:
		//static methods
		static std::map<size_t, std::string> static_meteoparamname; ///<Associate a name with meteo parameters in Parameters
		static std::vector<std::string> s_default_paramname;
		static const double epsilon; ///<for comparing fields
		static const bool __init;    ///<helper variable to enable the init of static collection data
		static bool initStaticData();///<initialize the static map meteoparamname

		//private data members, please keep the order consistent with declaration lists and logic!
		std::vector<std::string> param_name;
		std::vector<double> data;
		size_t nrOfAllParameters;
		bool resampled; ///<set this to true if MeteoData is result of resampling
};

} //end namespace

#endif
