/***********************************************************************************/
/*  Copyright 2012 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef __PROCEXPSMOOTHING_H__
#define __PROCEXPSMOOTHING_H__

#include <meteoio/meteofilters/WindowedFilter.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @brief Exponential smoothing.
 * This implements an exponential moving average smoothing of parameter alpha (https://en.wikipedia.org/wiki/Exponential_smoothing)
 * such as (s being the filtered values and x the raw data):
 * - s_0 = x_0
 * - s_n = alpha*x_(t-1) + (1-alpha)*s_(t-1)
 *
 * Nodata values are excluded from the moving average calculation. It takes the following options:
 * - an optional "soft" keyword to allow (or not) the window centering to be adjusted to the available data
 * - the keywords "left", "center" or "right", indicating the window centering strategy
 * - the minimal number of points in window
 * - the minimal time interval spanning the window (in seconds)
 * - the alpha parameter, such as 0 < alpha < 1
 *
 * The standard behavior for this filter is obtained by using a left window. If using a right window, it behaves as if time was reversed
 * (ie. predictions from the future). A centered window applies the standard algorithms on the <b>distance</b> between the center point and each points,
 * that is that points get averaged symetrically around the middle point before running the standard algorithm.
 *
 * @note This filter temporally shifts the signal back by window_width/2 if using left or right window
 * @note This would probably lead to slightly unexpected results if used on irregularly sampled data
 *
 * @code
 *          TA::filter1 = exp_smoothing
 *          TA::arg1    = right 1 1800 0.8 ;(1800 seconds time span for the strictly right leaning window), alpha=0.8
 *          RH::filter1 = mean_avg
 *          RH::arg1    = 10 600 0.6 ;(strictly left window spanning 600 seconds and at least 10 points), alpha=0.6
 * @endcode
 */

class ProcExpSmoothing : public WindowedFilter {
	public:
		ProcExpSmoothing(const std::vector<std::string>& vec_args, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(std::vector<std::string> vec_args);
		double calcExpSmoothing(const std::vector<MeteoData>& ivec, const unsigned int& param, const size_t& start, const size_t& end, const size_t& pos);

		double alpha;
};

} //end namespace

#endif
