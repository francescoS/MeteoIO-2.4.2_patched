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

#ifndef __METEO2DINTERPOLATOR_H__
#define __METEO2DINTERPOLATOR_H__

#include <meteoio/Config.h>
#include <meteoio/Date.h>
#include <meteoio/MeteoData.h>
#include <meteoio/DEMObject.h>
#include <meteoio/InterpolationAlgorithms.h>

#include <memory>
#include <vector>
#include <map>

namespace mio {

class IOManager; //forward declaration

/**
 * @page dev_2Dinterpol How to write a spatial interpolation algorithm
 * Point measurements can be spatially interpolated by MeteoIO, through the use of interpolation
 * algorithms. The user will then be able to choose for each meteorological parameter which
 * interpolations could be applicable and the system will choose (through a heuristic) which one
 * to apply at each time step (depending on the conditions of the moment, like the number of measurements).
 *
 * @section structure_2Dinterpol Structure
 * The selection of which interpolation algorithm to use at any given time step, for a given parameter is
 * performed by the Meteo2DInterpolator class. This class provides the interface for the spatial
 * interpolations. The interpolation algorithms themselves derive from the
 * InterpolationAlgorithm class that standardizes their public API (which would not be used by anything
 * but a Meteo2DInterpolator object). It contains a getQualityRating() method that must return a quality
 * index between 0 (algorithm not applicable) and 1 (perfect result if using this algorithm). This is currently
 * only based on extremely simple heuristics, using general knowledge about the applicability of the various
 * spatial interpolation methods depending on some obvious factors (number of measurement points, etc). The
 * Meteo2DInterpolator object will call this method from all the algorithms listed by the user (in his io.ini
 * configuration file) and keep the one that gets the highest score for interpolating the current parameter
 * at the current time step. The interpolation is then done calling the algorithm's calculate method.
 *
 * @section implementation_2Dinterpol Implementation
 * It is therefore necessary to create in InterpolationAlgorithms.cc (and declared in the .h) a new class,
 * nammed after the algorithm that will be implemented and inheriting InterpolationAlgorithm. Three methods need
 * to be implemented (the constructor being inherited from InterpolationAlgorithm and automatically called
 * by an object factory):
 * - void initialize(const MeteoData::Parameters& in_param)
 * - double getQualityRating()
 * - void calculate(Grid2DObject& grid)
 *
 * The initialize method takes the meteorological parameter that will be interpolated and set the param
 * private member to it. It then computes the private member nrOfMeasurments that contains the number of
 * stations that have this meteorological parameter available by calling getData(param, vecData, vecMeta), which
 * also fills the vectors vecData and vecMeta with the available data (as double) and metadata (as StationData).
 * Custom data preparation can obviously be done in this method.
 *
 * The calculate method must properly erase and reste the grid that it receives before filling it. If necessary,
 * (as is the case for precipitation, relative humidity and snow height, for example) the grid can be checked for min/max by
 * calling checkMinMax() at the end of Meteo2DInterpolator::interpolate.It can also add extra information about the
 * interpolation process (such as a regression coefficient or error estimate) to the InterpolationAlgorithm::info
 * stringstream (which will be made available to external programs, such as GUIs).
 *
 * The new class and its associated end user key must be used and its constructor called in AlgorithmFactory::getAlgorithm.
 * It is recommended that any generic statistical
 * spatial processing be implemented as a static class in libinterpol2D.cc so that it could be reused by other
 * algorithms (see for example Interpol2D::IDW and IDWCore). In any case, proper doxygen documentation
 * must be written alongside the implementation.
 *
 * @section doc_2Dinterpol Documentation
 * The newly added interpolation algorithm must be added to the list of available algorithms in
 * InterpolationAlgorithms.h with a proper description. An example can also be given in the example section
 * of the same file. Please feel free to add necessary bibliographic references to the bibliographic section!
 *
*/

/**
 * @class Meteo2DInterpolator
 * @brief A class to spatially interpolate meteo parameters. For more, see \ref interpol2d
 *
 * @ingroup stats
 * @author Mathias Bavay and Thomas Egger
 * @date   2010-01-14
 */

#ifdef _POPC_
#include <paroc_base.h>
class Meteo2DInterpolator : POPBase {
	public:
		void Serialize(POPBuffer &buf, bool pack);
#else
class Meteo2DInterpolator {
#endif
 	public:
		/**
		* @brief Constructor.
		*/
		Meteo2DInterpolator(const Config& i_cfg);
		Meteo2DInterpolator(const Config& i_cfg, IOManager& iomanager);
		Meteo2DInterpolator(const Meteo2DInterpolator& c);

		~Meteo2DInterpolator();

		/**
		 * @brief A generic function that can interpolate for any given MeteoData member variable
		 *
		 * @param date date for which to interpolate
		 * @param dem Digital Elevation Model on which to perform the interpolation
		 * @param meteoparam Any MeteoData member variable as specified in the
		 * 				 enum MeteoData::Parameters (e.g. MeteoData::TA)
		 * @param result A Grid2DObject that will be filled with the interpolated data
		 */
		void interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
		                 Grid2DObject& result);

		/**
		 * @brief A generic function that can interpolate for any given MeteoData member variable
		 *
		 * @param date date for which to interpolate
		 * @param dem Digital Elevation Model on which to perform the interpolation
		 * @param meteoparam Any MeteoData member variable as specified in the
		 * 				 enum MeteoData::Parameters (e.g. MeteoData::TA)
		 * @param result A Grid2DObject that will be filled with the interpolated data
		 * @param InfoString some information about the interpolation process (useful for GUIs)
		 */
		void interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
		                 Grid2DObject& result, std::string& InfoString);

		void interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                            const std::vector<Coords>& in_coords, const bool& use_full_dem, std::vector<double>& result, std::string& info_string);

		/**
		 * @brief Retrieve the arguments vector for a given interpolation algorithm
		 * @param param the meteorological parameter that is concerned
		 * @param algorithm the desired algorithm
		 * @param vecArgs a vector of strings containing the arguments
		 */
		size_t getArgumentsForAlgorithm(const std::string& param,
		                                const std::string& algorithm,
		                                std::vector<std::string>& vecArgs) const;

		void setIOManager(IOManager& iomanager);
		Meteo2DInterpolator& operator=(const Meteo2DInterpolator& source);
		const std::string toString() const;

	private:
		static void checkMinMax(const double& minval, const double& maxval, Grid2DObject& gridobj);
		static void check_projections(const DEMObject& dem, const std::vector<MeteoData>& vec_meteo);
		static size_t get_parameters(const Config& cfg, std::set<std::string>& set_parameters);
		static size_t getAlgorithmsForParameter(const Config& cfg, const std::string& parname, std::vector<std::string>& vecAlgorithms);

		void addToBuffer(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam, const Grid2DObject& grid);
		bool getFromBuffer(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam, Grid2DObject& grid) const;
		void setDfltBufferProperties();
		void setAlgorithms();

		const Config& cfg; ///< Reference to Config object, initialized during construction
		IOManager *iomanager; ///< Reference to IOManager object, used for callbacks, initialized during construction
		std::map<std::string, Grid2DObject> mapBufferedGrids; ///< Buffer interpolated grids
		std::vector<std::string> IndexBufferedGrids; ///< Keep position information for easy erase fo specific grids
		std::map< std::string, std::vector<InterpolationAlgorithm*> > mapAlgorithms; //per parameter interpolation algorithms
		size_t max_grids; ///< How many grids to buffer
		bool algorithms_ready; ///< Have the algorithms objects been constructed?
};

} //end namespace

#endif
