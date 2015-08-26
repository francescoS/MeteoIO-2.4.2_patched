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
#include "PGMIO.h"
#include <errno.h>
#include <string.h>
#include <algorithm>

using namespace std;

namespace mio {
/**
 * @page pgmio PGMIO
 * @section pgmio_format Format
 * This reads a grid file in PGM format (see http://www.fileformat.info/format/pbm/egff.htm). This is a graphic format that is supported by a wide range of graphics programs (Gimp, Irfanview, Paint Shop Pro, gqview, etc). This allows to write a grid as an image (one pixel equals one cell), read an image as a grid (useful for creating synthetic DEMs). Since there is no geolocalization information in this format, such data is either encoded as a comment (when writing a file) or read from io.ini (for reading).
 * Finally, the naming scheme for meteo grids should be: YYYY-MM-DDTHH.mm_{MeteoGrids::Parameters}.pgm
 *
 * Please keep in mind that only a finite number of greyscales are used, making a discretization of the data. Moreover, we consider that a color of "0" is NODATA.
 *
 * @section pgmio_units Units
 * Cellsize is in meters, x/y coords in the section's coordinate system.
 *
 * @section pgmio_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - GRID2DPATH: meteo grids directory where to read/write the grids; [Input] and [Output] sections
 * - PGM_XCOORD: lower left x coordinate; [Input] section
 * - PGM_YCOORD: lower left y coordinate; [Input] section
 * - PGM_CELLSIZE: cellsize in meters; [Input] section
 * - PGM_MIN: minimum value in real world coordinates to match with the minimum value read out of the PGM file (such minimum being greater than 0 because 0 is NODATA)
 * - PGM_MAX: maximum value in real world coordinates to match with the maximum value read out of the PGM file
 */

const double PGMIO::plugin_nodata = 0.; //plugin specific nodata value. It can also be read by the plugin (depending on what is appropriate)

PGMIO::PGMIO(const std::string& configfile)
       : cfg(configfile),
         coordin(), coordinparam(), coordout(), coordoutparam(),
         fin(), fout(), grid2dpath_in(), grid2dpath_out()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	getGridPaths();
}

PGMIO::PGMIO(const Config& cfgreader)
       : cfg(cfgreader),
         coordin(), coordinparam(), coordout(), coordoutparam(),
         fin(), fout(), grid2dpath_in(), grid2dpath_out()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	getGridPaths();
}

void PGMIO::getGridPaths() {
	grid2dpath_in.clear(), grid2dpath_out.clear();
	string tmp;
	cfg.getValue("GRID2D", "Input", tmp, IOUtils::nothrow);
	if (tmp == "PGM") //keep it synchronized with IOHandler.cc for plugin mapping!!
		cfg.getValue("GRID2DPATH", "Input", grid2dpath_in);
	tmp.clear();
	cfg.getValue("GRID2D", "Output", tmp, IOUtils::nothrow);
	if (tmp == "PGM") //keep it synchronized with IOHandler.cc for plugin mapping!!
		cfg.getValue("GRID2DPATH", "Output", grid2dpath_out);
}

PGMIO::~PGMIO() throw() {

}

size_t PGMIO::getNextHeader(std::vector<std::string>& vecString, const std::string& filename) {
	std::string line;

	while(!fin.eof()) {
		getline(fin, line);
		IOUtils::trim(line);
		if(!line.empty() && line.at(0)!='#') {
			return IOUtils::readLineToVec(line, vecString);
		}
	}
	throw IOException("Can not read necessary header lines in " + filename, AT);
}

void PGMIO::read2DGrid_internal(Grid2DObject& grid_out, const std::string& full_name)
{
	size_t ncols, nrows;
	unsigned int nr_colors;
	double xllcorner, yllcorner, cellsize;
	double tmp_val, val_min, val_max;
	std::vector<std::string> tmpvec;
	std::string line;

	if (!IOUtils::validFileName(full_name)) {
		throw InvalidFileNameException(full_name, AT);
	}
	if (!IOUtils::fileExists(full_name)) {
		throw FileNotFoundException(full_name, AT);
	}

	fin.clear();
	fin.open (full_name.c_str(), ifstream::in);
	if (fin.fail()) {
		ostringstream ss;
		ss << "Error opening file \"" << full_name << "\", possible reason: " << strerror(errno);
		throw FileAccessException(ss.str(), AT);
	}

	const char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

	//Go through file, save key value pairs
	try {
		//read header: magic value
		if(getNextHeader(tmpvec, full_name)!=1) {
			throw IOException("Can not read necessary header in " + full_name, AT);
		}
		//read header: image width and height
		if(getNextHeader(tmpvec, full_name)!=2) {
			throw IOException("Can not read necessary header in " + full_name, AT);
		}
		IOUtils::convertString(ncols, tmpvec[0]);
		IOUtils::convertString(nrows, tmpvec[1]);
		//read header: number of greys
		if(getNextHeader(tmpvec, full_name)!=1) {
			throw IOException("Can not read necessary header in " + full_name, AT);
		}
		IOUtils::convertString(nr_colors, tmpvec[0]);

		cfg.getValue("PGM_XCOORD", "Input", xllcorner, IOUtils::dothrow);
		cfg.getValue("PGM_YCOORD", "Input", yllcorner, IOUtils::dothrow);
		cfg.getValue("PGM_CELLSIZE", "Input", cellsize, IOUtils::dothrow);
		cfg.getValue("PGM_MIN", "Input", val_min, IOUtils::dothrow);
		cfg.getValue("PGM_MAX", "Input", val_max, IOUtils::dothrow);

		Coords location(coordin, coordinparam);
		location.setXY(xllcorner, yllcorner, IOUtils::nodata);

		//Initialize the 2D grid
		grid_out.set(ncols, nrows, cellsize, location);

		//initialize scale factor
		const double scale_factor = (val_max-val_min)/(double)(nr_colors-2); //because 256 colors = 0 to 255!! and color0 = nodata

		//Read one line after the other and parse values into Grid2DObject
		for (size_t kk=nrows-1; (kk < nrows); kk--) {
			getline(fin, line, eoln); //read complete line

			if (IOUtils::readLineToVec(line, tmpvec) != ncols) {
				ostringstream ss;
				ss << "Invalid number of columns at line " << nrows-kk << " in file \"" << full_name << "\". ";
				ss << "Expecting " << ncols << " columns\n";
				throw InvalidFormatException(ss.str(), AT);
			}

			for (size_t ll=0; ll < ncols; ll++){
				if (!IOUtils::convertString(tmp_val, tmpvec[ll], std::dec)) {
					throw ConversionFailedException("For Grid2D value in line: " + line + " in file " + full_name, AT);
				}

				if(tmp_val==plugin_nodata) {
					//replace file's nodata by uniform, internal nodata
					grid_out(ll, kk) = IOUtils::nodata;
				} else {
					grid_out(ll, kk) = (tmp_val-1)*scale_factor+val_min; //because color0 = nodata
				}
			}
		}
	} catch(const std::exception&) {
		cerr << "[E] error when reading PGM grid \"" << full_name << "\" " << AT << ": "<< endl;
		cleanup();
		throw;
	}
	cleanup();
}

void PGMIO::read2DGrid(Grid2DObject& grid_out, const std::string& filename) {
	read2DGrid_internal(grid_out, grid2dpath_in+"/"+filename);
}

void PGMIO::read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	std::string date_str = date.toString(Date::ISO);
	std::replace( date_str.begin(), date_str.end(), ':', '.');
	read2DGrid_internal(grid_out, grid2dpath_in+"/"+date_str+"_"+MeteoGrids::getParameterName(parameter)+".pgm");
}

void PGMIO::readDEM(DEMObject& dem_out)
{
	string filename;
	cfg.getValue("DEMFILE", "Input", filename);
	read2DGrid_internal(dem_out, filename);
}

void PGMIO::readLanduse(Grid2DObject& /*landuse_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PGMIO::readAssimilationData(const Date& /*date_in*/, Grid2DObject& /*da_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PGMIO::readStationData(const Date&, std::vector<StationData>& /*vecStation*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PGMIO::readMeteoData(const Date& /*dateStart*/, const Date& /*dateEnd*/,
                          std::vector< std::vector<MeteoData> >& /*vecMeteo*/,
                          const size_t&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PGMIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& /*vecMeteo*/,
                           const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PGMIO::readPOI(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void PGMIO::write2DGrid(const Grid2DObject& grid_in, const std::string& name)
{
	std::string full_name = grid2dpath_out+"/"+name;
	const unsigned int nr_colors = 256;
	fout.open(full_name.c_str());
	if (fout.fail()) {
		ostringstream ss;
		ss << "Error opening file \"" << full_name << "\", possible reason: " << strerror(errno);
		throw FileAccessException(ss.str(), AT);
	}

	Coords llcorner=grid_in.llcorner;
	//we want to make sure that we are using the provided projection parameters
	//so that we output is done in the same system as the inputs
	llcorner.setProj(coordout, coordoutparam);

	fout << fixed << showpoint << setprecision(6);

	try {
		const double min_value = grid_in.grid2D.getMin();
		const double max_value = grid_in.grid2D.getMax();
		double scaling;
		if(min_value!=max_value) scaling = 1./(max_value - min_value) * (double)(nr_colors-1); //so we keep color 0 for nodata
		else scaling = 1.;

		//writing the header
		fout << "P2\n";
		fout << "#Generated by MeteoIO - http://slfsmm.indefero.net/p/meteoio\n";
		fout << "#llcorner latitude = " << setprecision(6) << llcorner.getLat() << "\n";
		fout << "#llcorner longitude = " << setprecision(6) << llcorner.getLon() << "\n";
		fout << "#cellsize = " << setprecision(2) << grid_in.cellsize << " m\n";
		fout << "#minimum = " << setprecision(6) << min_value << "\n";
		fout << "#maximum = " << setprecision(6) << max_value << "\n";
		fout << grid_in.ncols << " " << grid_in.nrows << "\n";
		fout << nr_colors << "\n";

		//writing the data
		if(grid_in.nrows>0) {
			for (size_t kk=grid_in.nrows-1; kk < grid_in.nrows; kk--) {
				for (size_t ll=0; ll < grid_in.ncols; ll++) {
					const double value = grid_in(ll, kk);
					if(value!=IOUtils::nodata)
						fout << static_cast<unsigned int>( floor((grid_in(ll, kk)-min_value)*scaling)+1 ) << " ";
					else
						fout << "0" << " ";
				}
				fout << "\n";
			}
		}
	} catch(...) {
		cerr << "[E] error when writing PGM grid \"" << full_name << "\" " << AT << ": "<< endl;
		cleanup();
		throw;
	}

	cleanup();
}

void PGMIO::write2DGrid(const Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	//path will be added by calling write2Dgrid
	std::string date_str = date.toString(Date::ISO);
	std::replace( date_str.begin(), date_str.end(), ':', '.');
	write2DGrid(grid_out, date_str+"_"+MeteoGrids::getParameterName(parameter)+".pgm");
}

void PGMIO::cleanup() throw()
{
	if (fin.is_open()) {//close fin if open
		fin.close();
	}
	if (fout.is_open()) {//close fout if open
		fout.close();
	}
}

} //namespace
