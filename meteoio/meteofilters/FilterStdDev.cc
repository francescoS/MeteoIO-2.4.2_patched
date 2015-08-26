/***********************************************************************************/
/*  Copyright 2011 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/meteofilters/FilterStdDev.h>
#include <cmath>

using namespace std;

namespace mio {

const double FilterStdDev::sigma = 2.; ///<How many times the stddev allowed for valid points

FilterStdDev::FilterStdDev(const std::vector<std::string>& vec_args, const std::string& name) : WindowedFilter(name)
{
	parse_args(vec_args);

	//This is safe, but maybe too imprecise:
	properties.time_before = min_time_span;
	properties.time_after  = min_time_span;
	properties.points_before = min_data_points;
	properties.points_after = min_data_points;
}

void FilterStdDev::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                           std::vector<MeteoData>& ovec)
{
	ovec = ivec;
	for (size_t ii=0; ii<ovec.size(); ii++){ //for every element in ivec, get a window
		double& value = ovec[ii](param);
		if(value==IOUtils::nodata) continue;

		//Calculate deviation
		double mean     = IOUtils::nodata;
		double std_dev  = IOUtils::nodata;

		size_t start, end;
		if( get_window_specs(ii, ivec, start, end) ) {
			getStat(ivec, param, start, end, std_dev, mean);
		}  else if(!is_soft) {
			value = IOUtils::nodata;
			continue;
		}

		if(mean==IOUtils::nodata) {
			if(!is_soft) value = IOUtils::nodata;
			continue;
		}

		if( value!=IOUtils::nodata && abs(value-mean)>sigma*std_dev) {
			value = IOUtils::nodata;
		}
	}
}

void FilterStdDev::getStat(const std::vector<MeteoData>& ivec, const unsigned int& param,
                           const size_t& start, const size_t& end, double& stddev, double& mean)
{
	size_t count=0;
	double sum=0.;

	for(size_t ii=start; ii<=end; ii++) {
		const double& value = ivec[ii](param);
		if(value!=IOUtils::nodata) {
			sum += value;
			count++;
		}
	}

	if(count<=1) {
		mean = IOUtils::nodata;
		stddev = IOUtils::nodata;
	} else {
		//compensated variance algorithm, see https://secure.wikimedia.org/wikipedia/en/wiki/Algorithms_for_calculating_variance
		mean = sum/(double)count;
		double sum2=0., sum3=0.;
		for(size_t ii=start; ii<=end; ii++) {
			const double& value = ivec[ii](param);
			if(value!=IOUtils::nodata) {
				const double delta = value - mean;
				sum2 += delta*delta;
				sum3 += delta;
			}
		}
		const double variance = (sum2 - sum3*sum3/static_cast<double>(count)) / static_cast<double>(count - 1);
		stddev = sqrt(variance);
	}
}

void FilterStdDev::parse_args(std::vector<std::string> vec_args)
{
	vector<double> filter_args;

	if (vec_args.size() > 2){
		is_soft = ProcessingBlock::is_soft(vec_args);
	}

	if (vec_args.size() > 2)
		centering = (WindowedFilter::Centering)WindowedFilter::get_centering(vec_args);

	convert_args(2, 2, vec_args, filter_args);

	if ((filter_args[0] < 1) || (filter_args[1] < 0)){
		throw InvalidArgumentException("Invalid window size configuration for filter " + getName(), AT);
	}

	min_data_points = (unsigned int)floor(filter_args[0]);
	min_time_span = Duration(filter_args[1] / 86400.0, 0.);
}

}
