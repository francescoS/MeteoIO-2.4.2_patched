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
#ifndef __FILTERMEANAVG_H__
#define __FILTERMEANAVG_H__

#include <meteoio/meteofilters/WindowedFilter.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  FilterMeanAvg
 * @ingroup processing
 * @author Thomas Egger
 * @date   2011-01-24
 * @brief Mean average processing.
 * The mean average filter returns the mean value of all values within a user given time window. Remarks:
 * - nodata values are excluded from the mean
 * - Two arguments expected (both have to be fullfilled for the filter to start operating):
 *   - minimal number of points in window
 *   - minimal time interval spanning the window (in seconds)
 * - the two arguments may be preceded by the keywords "left", "center" or "right", indicating the window position
 * - the keyword "soft" maybe added, if the window position is allowed to be adjusted to the data present
 *
 * @code
 * Valid examples for the io.ini file:
 *          TA::filter1 = mean_avg
 *          TA::arg1    = soft left 1 1800 (1800 seconds time span for the left leaning window)
 *          RH::filter1 = mean_avg
 *          RH::arg1    = 10 600          (strictly centered window spanning 600 seconds and at least 10 points)
 * @endcode
 */

class FilterMeanAvg : public WindowedFilter {
	public:
		FilterMeanAvg(const std::vector<std::string>& vec_args, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(std::vector<std::string> vec_args);
		double calc_avg(const std::vector<MeteoData>& ivec, const unsigned int& param, const size_t& start, const size_t& end);
};

} //end namespace

#endif
