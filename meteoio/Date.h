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
#ifndef __DATE_H__
#define __DATE_H__

#include <meteoio/IOExceptions.h>

#include <string>
#include <sstream>

//#define NEGATIVE_JULIAN
namespace mio {

/**
 * @class Date
 * @brief  A class to handle timestamps.
 * This class handles conversion between different time display formats (ISO, numeric) as well as different
 * time representation (julian date, modified julian date, etc). It also handles time zones as well as
 * very basic Daylight Saving Time (DST). Since the activation dates of DST are political and not technical,
 * it can not be automatically calculated. Therefore, it has to be provided by the caller: when the dst flag
 * is set, the dst time shift is automatically applied. When the dst flag ceases to be set, the dst time shift
 * is no longer applied. This is very crude, but please keep in mind that using DST for monitoring data is
 * usually a bad idea... Finally, we assume that dates are positive. If this would not be the case, this
 * class has to be recompiled with the proper define.
 *
 * Internally, the date is stored as true julian date in GMT.
 * The maximal precision is 1 minute (that can be easily brought to 1 seconds if
 * it would appear necessary/useful, with the limitation that leap seconds are currently not handled).
 *
 * Please see Date::FORMATS for supported display formats and http://en.wikipedia.org/wiki/Julian_day for
 * the various date representation definitions. The following data representation are currently supported:
 * - julian date, see Date::getJulianDate
 * - modified julian date, see Date::getModifiedJulianDate
 * - truncated julian date, see Date::getTruncatedJulianDate
 * - Unix date, see Date::getUnixDate
 * - Excel date, see Date::getExcelDate
 *
 * @ingroup data_str
 * @author Mathias Bavay
 * @date 2010-04-15
 */

#ifdef _POPC_
#include <paroc_base.h>
class DateDummy {}; //HACK for POPC

class Date : POPBase {
	public:
		void Serialize(POPBuffer &buf, bool pack);
#else
class Date {
#endif
	public:
		///Keywords for selecting the date formats
		typedef enum {
			ISO, ///< ISO 8601 extended format combined date: YYYY-MM-DDTHH:mm:SS (fields might be dropped, in the least to the most significant order)
			ISO_TZ, ///< ISO 8601 format (same as ISO) but with time zone specification
			FULL, ///< ISO 8601 followed by the julian date (in parenthesis)
			NUM, ///< ISO 8601 basic format date: YYYYMMDDHHmmSS (fields might be dropped, in the least to the most significant order)
			DIN ///<DIN5008 format: DD.MM.YYYY HH:MM
		} FORMATS;

		///Keywords for selecting rounding strategy
		typedef enum RND_TYPE {
			UP, ///< rounding toward highest absolute value
			DOWN, ///< rounding toward smallest absolute value
			CLOSEST ///< rounding toward closest
		} RND;

		static const int daysLeapYear[];
		static const int daysNonLeapYear[];
		static const double DST_shift;
		static const float MJD_offset;
		static const float Unix_offset;
		static const float Excel_offset;
		static const float Matlab_offset;

		Date();
		Date(const double& julian_in, const double& in_timezone, const bool& in_dst=false);
		Date(const int& year, const int& month, const int& day, const int& hour, const int& minute, const double& in_timezone, const bool& in_dst=false);
		Date(const time_t&, const bool& in_dst=false);

		void setFromSys();
		void setTimeZone(const double& in_timezone, const bool& in_dst=false);
		void setDate(const Date& in_date);
		void setDate(const double& julian_in, const double& in_timezone, const bool& in_dst=false);
		void setDate(const int& year, const int& month, const int& day, const int& hour, const int& minute, const double& in_timezone, const bool& in_dst=false);
		void setDate(const int& year, const unsigned int& month, const unsigned int& day, const unsigned int& hour, const unsigned int& minute, const double& in_timezone, const bool& in_dst=false);
		void setDate(const time_t& in_time, const bool& in_dst=false);
		void setModifiedJulianDate(const double& julian_in, const double& in_timezone, const bool& in_dst=false);
		void setUnixDate(const time_t& in_time, const bool& in_dst=false);
		void setExcelDate(const double excel_in, const double& in_timezone, const bool& in_dst=false);
		void setMatlabDate(const double excel_in, const double& in_timezone, const bool& in_dst=false);
		void setUndef(const bool& flag);

		bool isUndef() const;
		double getTimeZone() const;
		bool getDST() const;
		double getJulian(const bool& gmt=false) const;
		double getModifiedJulianDate(const bool& gmt=false) const;
		double getTruncatedJulianDate(const bool& gmt=false) const;
		time_t getUnixDate() const;
		double getExcelDate(const bool& gmt=false) const;
		double getMatlabDate(const bool& gmt=false) const;

		void getDate(double& julian_out, const bool& gmt=false) const;
		void getDate(int& year, int& month, int& day, const bool& gmt=false) const;
		void getDate(int& year, int& month, int& day, int& hour, const bool& gmt=false) const;
		void getDate(int& year, int& month, int& day, int& hour, int& minute, const bool& gmt=false) const;
		int getYear(const bool& gmt=false) const;

		int getJulianDayNumber(const bool& gmt=false) const;
		bool isLeapYear() const;

		static double rnd(const double& julian, const unsigned int& precision, const RND& type=CLOSEST);
		void rnd(const unsigned int& precision, const RND& type=CLOSEST);
		static const Date rnd(const Date& indate, const unsigned int& precision, const RND& type=CLOSEST);
		static double parseTimeZone(const std::string& timezone_iso);

		static std::string printFractionalDay(const double& fractional);
		const std::string toString(FORMATS type, const bool& gmt=false) const;
		const std::string toString() const;
		friend std::iostream& operator<<(std::iostream& os, const Date& date);
		friend std::iostream& operator>>(std::iostream& is, Date& date);

		//Operator Prototypes
		bool operator==(const Date&) const;
		bool operator!=(const Date&) const;
		bool operator<(const Date&) const;
		bool operator<=(const Date&) const;
		bool operator>(const Date&) const;
		bool operator>=(const Date&) const;

		///Intervals arithmetic
		///Can be used to add an interval to an existing Date object.
		///Construct a Date object representing the interval e.g. Date(1.0) for 1 day and add that to another Date object.
		///Please use the Duration type instead of Date for such calculations!
		Date& operator+=(const Date&);
		Date& operator-=(const Date&);
		Date& operator+=(const double&);
		Date& operator-=(const double&);
		Date& operator*=(const double&);
		Date& operator/=(const double&);

		const Date operator+(const Date&) const;
		const Date operator-(const Date&) const;
		const Date operator+(const double&) const;
		const Date operator-(const double&) const;
		const Date operator*(const double&) const;
		const Date operator/(const double&) const;

	protected:
		double localToGMT(const double& in_julian) const;
		double GMTToLocal(const double& in_gmt_julian) const;
		double calculateJulianDate(const int& in_year, const int& in_month, const int& in_day, const int& in_hour, const int& in_minute) const;
		void calculateValues(const double& i_julian, int& out_year, int& out_month, int& out_day, int& out_hour, int& out_minute) const;
		long getJulianDayNumber(const int&, const int&, const int&) const;
		bool isLeapYear(const int&) const;
		void plausibilityCheck(const int& in_year, const int& in_month, const int& in_day, const int& in_hour, const int& in_minute) const;

		static const double epsilon;
		double timezone;
		double gmt_julian;
		int gmt_year, gmt_month, gmt_day, gmt_hour, gmt_minute;
		bool dst;
		bool undef;
};

typedef Date Duration; //so that later, we can implement a true Interval/Duration class

} //end namespace

#endif
