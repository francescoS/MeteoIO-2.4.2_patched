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
#ifndef __IOUTILS_H__
#define __IOUTILS_H__

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <cstdlib>
#include <limits>
#include <cmath>

#include <meteoio/Coords.h>
#include <meteoio/Date.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/meteolaws/Meteoconst.h>

#ifndef C_TO_K
#define C_TO_K( T ) ( T + Cst::t_water_freezing_pt )	  // degree Celsius to kelvin
#endif
#ifndef K_TO_C
#define K_TO_C( T ) ( T - Cst::t_water_freezing_pt )	  // kelvin to degree Celsius
#endif

namespace mio {

#ifdef _MSC_VER
double round(const double& x);
#endif

class MeteoData;
class Coords;
class Config;

/**
* @brief Return the library version
* @return library version string
*/
std::string getLibVersion();

namespace IOUtils {
	enum ThrowOptions { dothrow, nothrow };
	const double nodata = -999.0; ///<This is the internal nodata value
	//const double not_set = std::numeric_limits<double>::max()-2.;
	const unsigned int unodata = (unsigned int)-1;
	const int inodata = -999;
	const short int snodata = -999;
	const size_t npos    = (size_t)-1;  ///<npos is the out-of-range value

	const double grid_epsilon = 5.; ///<What is an acceptable small distance on a grid, in meters
	const double lon_epsilon = grid_epsilon / Cst::earth_R0; ///<in degrees. Small angle for longitudes, so sin(x)=x
	const double lat_epsilon = lon_epsilon/2.; ///<in degrees. Small angle for latitudes. Since for longitudes it is for 360deg, it has to be 180deg for latitudes

	/**
	* @brief Check whether two values are equal regarding a certain epsilon environment (within certain radius of each other)
	* @param val1
	* @param val2
	* @param epsilon is a radius around val1
	* @return true if val2 is within the radius around val1, false otherwise.
	*/
	inline bool checkEpsilonEquality(const double& val1, const double& val2, const double& epsilon) {return (fabs(val1-val2) < epsilon);}

	size_t seek(const Date& soughtdate, const std::vector<MeteoData>& vecM, const bool& exactmatch=true);

	/**
	* @brief Converts a compass bearing to a trigonometric angle
	* @param bearing compass bearing (0° on top, clockwise, in [0°, 360°[)
	* @return trigonometric angle (0° on the right, counterclockwise, in [0, 2PI[)
	*/
	double bearing_to_angle(const double& bearing);
	/**
	* @brief Converts a trigonometric angle to a compass bearing
	* @param angle trigonometric angle (0° on the right, counterclockwise, in [0, 2PI[)
	* @return bearing (0° on top, clockwise, in [0°, 360°[)
	*/
	double angle_to_bearing(const double& angle);
	/**
	* @brief Converts a string bearing to a compass bearing
	* @param bearing_str as N, NE, SSW, etc
	* @return bearing (0° on top, clockwise, in [0°, 360°[)
	*/
	double bearing(std::string bearing_str);

	/**
	* @brief Retrieve the user name
	* This checks various environment variables (USERNAME, USER, LOGNAME).
	* @return user name
	*/
	std::string getLogName();

	/**
	* @brief Removes trailing and leading whitespaces, tabs and newlines from a string.
	* @param s The reference of the string to trim (in/out parameter)
	*/
	void trim(std::string &s);

	void stripComments(std::string& str);

	/**
	* @brief read a string line, parse it and save it into a map object, that is passed by reference
	* @param in_line (const string&) string to parse
	* @param delimiter (const string&) delimiter to use for the parsing
	* @param setToUpperCase If set to true the key will be put into upper case (for case insensitivity)
	* @param key retrieved key
	* @param value retrieved value
	* @return (bool) true when line is empty
	*/
	bool readKeyValuePair(const std::string& in_line, const std::string& delimiter,
	                      std::string &key, std::string &value, const bool& setToUpperCase=false);

	void toUpper(std::string& str);
	std::string strToUpper(std::string str);
	void toLower(std::string& str);
	std::string strToLower(std::string str);
	bool isNumeric(std::string input, const unsigned int& nBase=10);
	size_t readLineToVec(const std::string& line_in, std::vector<double>& vec_data);
	size_t readLineToVec(const std::string& line_in, std::vector<std::string>& vecString);
	size_t readLineToVec(const std::string& line_in, std::vector<std::string>& vecString, const char& delim);

	/**
	* @brief Convert a string to the requested type (template function).
	* @tparam T   [in] The type wanted for the return value (template type parameter).
	* @param t   [out] The value converted to the requested type.
	* @param str   [in] The input string to convert; trailling whitespaces are ignored,
	*              comment after non-string values are allowed, but multiple values are not allowed.
	* @param f   [in] The radix for reading numbers, such as std::dec or std::oct; default is std::dec.
	* @return true if everything went fine, false otherwise
	*/
	template <class T> bool convertString(T& t, const std::string& str, std::ios_base& (*f)(std::ios_base&) = std::dec) {
		std::string s(str);
		trim(s); //delete trailing and leading whitespaces and tabs
		if (s.empty()) {
			t = static_cast<T> (nodata);
			return true;
		} else {
			std::istringstream iss(s);
			iss.setf(std::ios::fixed);
			iss.precision(std::numeric_limits<double>::digits10); //try to read values with maximum precision
			iss >> f >> t; //Convert first part of stream with the formatter (e.g. std::dec, std::oct)
			if (iss.fail()) {
				//Conversion failed
				return false;
			}
			std::string tmp;
			getline(iss,  tmp); //get rest of line, if any
			trim(tmp);
			if (!tmp.empty() && tmp[0] != '#' && tmp[0] != ';') {
				//if line holds more than one value it's invalid
				return false;
			}
			return true;
		}
	}
	// fully specialized template functions (implementation must not be in header)
	template<> bool convertString<double>(double& t, const std::string& str, std::ios_base& (*f)(std::ios_base&));
	template<> bool convertString<std::string>(std::string& t, const std::string& str, std::ios_base& (*f)(std::ios_base&));
	template<> bool convertString<bool>(bool& t, const std::string& str, std::ios_base& (*f)(std::ios_base&));
	template<> bool convertString<unsigned int>(unsigned int& t, const std::string& str, std::ios_base& (*f)(std::ios_base&));
	template<> bool convertString<Coords>(Coords& t, const std::string& str, std::ios_base& (*f)(std::ios_base&));

	bool convertString(Date& t, const std::string& str, const double& time_zone, std::ios_base& (*f)(std::ios_base&) = std::dec);

	/**
	* @brief Returns, with the requested type, the value associated to a key (template function).
	* @tparam T   [in] The type wanted for the return value (template type parameter).
	* @param[in]  properties  A map containing all the parameters.
	* @param[in]  key         The key of the parameter to retrieve.
	* @param[out] t           The value associated to the key, converted to the requested type
	* @param[in]  options     Extra options, by default IOUtils::dothrow
	*/
	template <class T> void getValueForKey(const std::map<std::string,std::string>& properties,
	                                       const std::string& key, T& t, const ThrowOptions& options=IOUtils::dothrow){
		if (key.empty() && options!=IOUtils::nothrow)
			throw InvalidArgumentException("Empty key", AT);

		const std::map<std::string, std::string>::const_iterator it = properties.find(key);

		if (it == properties.end()){
			if (options == IOUtils::nothrow)
				return;
			else
				throw UnknownValueException("No value for key " + key, AT);
		}
		const std::string& value = it->second;

		if(!convertString<T>(t, value, std::dec) && options!=IOUtils::nothrow) {
			std::cerr << "[E] When reading \"" << key << "\" = \"" << t << "\"\n";
			throw ConversionFailedException(value, AT);
		}
	}

	/**
	* @brief Returns, with the requested type, the value associated to a key (template function).
	* @tparam T           [in] The type wanted for the return value (template type parameter).
	* @param[in]  properties  A map containing all the parameters.
	* @param[in]  key         The key of the parameter to retrieve.
	* @param[out] vecT        The vector of values associated to the key, each value is converted to the requested type
	* @param[in]  options     Extra options, by default IOUtils::dothrow
	*/
	template <class T> void getValueForKey(const std::map<std::string,std::string>& properties,
	                                       const std::string& key, std::vector<T>& vecT, const ThrowOptions& options=IOUtils::dothrow)
	{
		if (key.empty() && options!=IOUtils::nothrow)
			throw InvalidArgumentException("Empty key", AT);

		const std::map<std::string, std::string>::const_iterator it = properties.find(key);
		if (it == properties.end()) {
			if (options == IOUtils::nothrow) {
				return;
			} else {
				throw UnknownValueException("No value for key " + key, AT);
			}
		}
		const std::string& value = it->second;

		//split value string
		std::vector<std::string> vecUnconvertedValues;
		const size_t counter = readLineToVec(value, vecUnconvertedValues);
		for (size_t ii=0; ii<counter; ii++){
			T myvar;
			if(!convertString<T>(myvar, vecUnconvertedValues.at(ii), std::dec) && options!=IOUtils::nothrow){
				std::cerr << "[E] When reading \"" << key << "\" = \"" << myvar << "\"\n";
				throw ConversionFailedException(vecUnconvertedValues.at(ii), AT);
			}
			vecT.push_back(myvar);
		}
	}

	/**
	* @brief Standardize a given value to use MeteoIO's internal nodata value (if applicable)
	* @tparam T[in] The type wanted for the return value (template type parameter).
	* @param[in] value  the value to check/convert
	* @param[in] plugin_nodata plugin-specific nodata value
	* @return checked/converted value
	*/
	template <class T> T standardizeNodata(const T& value, const double& plugin_nodata) {
		if(value==plugin_nodata) return static_cast<T> (nodata);
		else return value;
	}

	/**
	* @brief A function that parses a Config object for COORSYS, COORDPARAM keywords in [Input] and [Output]
	*        section and sets the respective strings to the values of those keywords
	* @param[in] cfg  A Config object
	* @param[out] coordin The coordinate system to be used for input data
	* @param[out] coordinparam The coordinate system parameters to be used for output data
	* @param[out] coordout The coordinate system to be used for output data
	* @param[out] coordoutparam The coordinate system parameters to be used for output data
	*/
	void getProjectionParameters(const Config& cfg, std::string& coordin, std::string& coordinparam,
	                             std::string& coordout, std::string& coordoutparam);

	/**
	* @brief A function that parses a Config object for the time_zone keyword and returns the timezone
	* @param[in] cfg  A Config object
	* @param[out] tz_in value to be used for the input timezone
	* @param[out] tz_out value to be used for the output timezone
	*/
	void getTimeZoneParameters(const Config& cfg, double& tz_in, double& tz_out);

	/**
	* @brief Returns the parameters for splitting an array in several, balanced sub-arrays.
	* This is mostly usefull for parallel calculations, where an array will be split and sent to different
	* workers.
	* @param[in] dimx number of cells in the desired dimension
	* @param[in] nbworkers total number of slices
	* @param[in] wk current slice index (starting at 1)
	* @param[out] startx calculated start index for the current slice
	* @param[out] nx calculated number of cells (in the desired dimension) of the current slice
	*/
	void getArraySliceParams(const size_t& dimx, const unsigned int& nbworkers, const unsigned int &wk, size_t& startx, size_t& nx);

} //end namespace IOUtils

} //end namespace mio

#endif
