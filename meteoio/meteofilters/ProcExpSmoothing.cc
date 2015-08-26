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
#include <meteoio/meteofilters/ProcExpSmoothing.h>
#include <cmath>

using namespace std;

namespace mio {

ProcExpSmoothing::ProcExpSmoothing(const std::vector<std::string>& vec_args, const std::string& name)
                 : WindowedFilter(name), alpha(.5)
{
	parse_args(vec_args);

	//This is safe, but maybe too imprecise:
	properties.time_before = min_time_span;
	properties.time_after  = min_time_span;
	properties.points_before = min_data_points;
	properties.points_after = min_data_points;
}

void ProcExpSmoothing::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                            std::vector<MeteoData>& ovec)
{
	ovec = ivec;
	for (size_t ii=0; ii<ovec.size(); ii++){ //for every element in ivec, get a window
		double& value = ovec[ii](param);

		size_t start, end;
		if( get_window_specs(ii, ivec, start, end) ) {
			value = calcExpSmoothing(ivec, param, start, end, ii);
		} else if(!is_soft) value = IOUtils::nodata;
	}

}

double ProcExpSmoothing::calcExpSmoothing(const std::vector<MeteoData>& ivec, const unsigned int& param, const size_t& start, const size_t& end, const size_t& pos)
{
	const size_t max_len = max(pos-start, end-pos);
	bool initCompleted = false;
	double expavg = IOUtils::nodata;

	for (size_t ii=1; ii<=max_len; ii++) {
		//getting values left and right of the current point
		const double val1 = ( pos>=ii && (pos-ii)>=start )? ivec[pos-ii](param) : IOUtils::nodata;
		const double val2 = ( (pos+ii)<=end )? ivec[pos+ii](param) : IOUtils::nodata;

		//computing the average (centered window) or take the proper point (left or right window)
		double val;
		if(val1!=IOUtils::nodata && val2!=IOUtils::nodata) {
			val = (val1+val2) * .5;
		} else if(val1!=IOUtils::nodata) {
			val = val1;
		} else if(val2!=IOUtils::nodata) {
			val = val2;
		} else
			val = IOUtils::nodata;

		//add the contribution
		if(val!=IOUtils::nodata) {
			if (initCompleted){
				expavg = alpha*val + (1.-alpha)*expavg;
			} else {
				expavg = val;
				initCompleted = true;
			}
		}
	}

	return expavg;
}

void ProcExpSmoothing::parse_args(std::vector<std::string> vec_args)
{
	vector<double> filter_args;

	if (vec_args.size() > 2){
		is_soft = ProcessingBlock::is_soft(vec_args);
	}

	if (vec_args.size() > 2)
		centering = (WindowedFilter::Centering)WindowedFilter::get_centering(vec_args);

	convert_args(3, 3, vec_args, filter_args);

	if ((filter_args[0] < 1) || (filter_args[1] < 0)){
		throw InvalidArgumentException("Invalid window size configuration for filter " + getName(), AT);
	}

	min_data_points = (unsigned int)floor(filter_args[0]);
	min_time_span = Duration(filter_args[1] / 86400.0, 0.);
	alpha = filter_args[2];
	if(alpha<0. || alpha>1.) {
		throw InvalidArgumentException("The alpha parameter for filter " + getName() + " must be between 0 and 1!", AT);
	}
}

} //namespace
