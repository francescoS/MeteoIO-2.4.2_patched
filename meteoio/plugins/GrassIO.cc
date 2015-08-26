/***********************************************************************************/
/*  Copyright 2009 EPFL                                                            */
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
#include "GrassIO.h"

using namespace std;

namespace mio {
/**
 * @page grass GRASS
 * @section grass_format Format
 * This is for reading grid data in the JGRASS GIS format(see http://jgrass.wiki.software.bz.it/jgrass/JGrass_Wiki)
 *
 * @section grass_units Units
 * The distances are assumed to be in meters.
 *
 * @section grass_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: input coordinate system (see Coords) specified in the [Input] section
 * - COORDPARAM: extra input coordinates parameters (see Coords) specified in the [Input] section
 * - COORDSYS: output coordinate system (see Coords) specified in the [Output] section
 * - COORDPARAM: extra output coordinates parameters (see Coords) specified in the [Output] section
 * - DEMFILE: for reading the data as a DEMObject
 * - LANDUSE: for interpreting the data as landuse codes
 * - DAPATH: path+prefix of file containing data assimilation grids (named with ISO 8601 basic date and .sca extension, example ./input/dagrids/sdp_200812011530.sca)
 */

const double GrassIO::plugin_nodata = -999.0; //plugin specific nodata value

GrassIO::GrassIO(const std::string& configfile)
        : cfg(configfile), fin(), fout(), coordin(), coordinparam(), coordout(), coordoutparam()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

GrassIO::GrassIO(const Config& cfgreader)
        : cfg(cfgreader), fin(), fout(), coordin(), coordinparam(), coordout(), coordoutparam()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
}

GrassIO::~GrassIO() throw()
{
	cleanup();
}

void GrassIO::cleanup() throw()
{
	if (fin.is_open()) {//close fin if open
		fin.close();
	}
	if (fout.is_open()) {//close fout if open
		fout.close();
	}
}

void GrassIO::read2DGrid(Grid2DObject& grid_out, const std::string& filename)
{

	int _nx, _ny;
	size_t ncols, nrows;
	double north, east, south, west;
	double tmp_val, xllcorner, yllcorner, cellsize;
	vector<string> tmpvec;
	std::string line;
	std::map<std::string, std::string> header; // A map to save key value pairs of the file header

	if (!IOUtils::validFileName(filename)) {
		throw InvalidFileNameException(filename, AT);
	}
	if (!IOUtils::fileExists(filename)) {
		throw FileNotFoundException(filename, AT);
	}

	fin.clear();
	fin.open (filename.c_str(), ifstream::in);
	if (fin.fail()) {
		throw FileAccessException(filename, AT);
	}

	char eoln = IOUtils::getEoln(fin); //get the end of line character for the file

	//Go through file, save key value pairs
	try {
		IOUtils::readKeyValueHeader(header, fin, 6, ":");
		IOUtils::getValueForKey(header, "cols",  _nx);
		IOUtils::getValueForKey(header, "rows",  _ny);
		IOUtils::getValueForKey(header, "north", north);
		IOUtils::getValueForKey(header, "east",  east);
		IOUtils::getValueForKey(header, "south", south);
		IOUtils::getValueForKey(header, "west",  west);

		_nx = IOUtils::standardizeNodata(_nx, plugin_nodata);
		_ny = IOUtils::standardizeNodata(_ny, plugin_nodata);
		north = IOUtils::standardizeNodata(north, plugin_nodata);
		east = IOUtils::standardizeNodata(east, plugin_nodata);
		south = IOUtils::standardizeNodata(south, plugin_nodata);
		west = IOUtils::standardizeNodata(west, plugin_nodata);

		if ((_nx==0) || (_ny==0)) {
			throw IOException("Number of rows or columns in 2D Grid given is zero, in file: " + filename, AT);
		}
		if((_nx<0) || (_ny<0)) {
			throw IOException("Number of rows or columns in 2D Grid read as \"nodata\", in file: " + filename, AT);
		}
		ncols = (size_t)_nx;
		nrows = (size_t)_ny;
		xllcorner = west;
		yllcorner = south;
		cellsize = (east - west) / (double)ncols;

		//compute WGS coordinates (considered as the true reference)
		Coords coordinate(coordin, coordinparam);
		coordinate.setXY(xllcorner, yllcorner, IOUtils::nodata);

		//Initialize the 2D grid
		grid_out.set(ncols, nrows, cellsize, coordinate);

		//Read one line after the other and parse values into Grid2DObject
		for (size_t kk=nrows-1; (kk < nrows); kk--) {
			getline(fin, line, eoln); //read complete line

			if (IOUtils::readLineToVec(line, tmpvec) != ncols) {
				throw InvalidFormatException("Premature End " + filename, AT);
			}

			for (size_t ll=0; ll < ncols; ll++){
				if (tmpvec[ll] == "*"){
					tmp_val = plugin_nodata;
				} else {
					if (!IOUtils::convertString(tmp_val, tmpvec[ll], std::dec)) {
						throw ConversionFailedException("For Grid2D value in line: " + line + " in file " + filename, AT);
					}
				}

				if(tmp_val <= plugin_nodata) {
					//replace file's nodata by uniform, internal nodata
					grid_out(ll, kk) = IOUtils::nodata;
				} else {
					grid_out(ll, kk) = tmp_val;
				}
			}
		}
	} catch(const std::exception&) {
		cleanup();
		throw;
	}
	cleanup();
}

void GrassIO::read2DGrid(Grid2DObject&, const MeteoGrids::Parameters&, const Date&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}


void GrassIO::readDEM(DEMObject& dem_out)
{
	string filename;
	cfg.getValue("DEMFILE", "Input", filename);
	read2DGrid(dem_out, filename);
}

void GrassIO::readLanduse(Grid2DObject& landuse_out)
{
	string filename;
	cfg.getValue("LANDUSEFILE", "Input", filename);
	read2DGrid(landuse_out, filename);
}

void GrassIO::readAssimilationData(const Date& date_in, Grid2DObject& da_out)
{
	int yyyy, MM, dd, hh, mm;
	date_in.getDate(yyyy, MM, dd, hh, mm);
	string filepath;

	cfg.getValue("DAPATH", "Input", filepath);

	ostringstream ss;
	ss.fill('0');
	ss << filepath << "/" << setw(4) << yyyy << setw(2) << MM << setw(2) <<  dd << setw(2) <<  hh << setw(2) <<  mm <<".sca";

	read2DGrid(da_out, ss.str());
}

void GrassIO::readStationData(const Date&, std::vector<StationData>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GrassIO::readMeteoData(const Date&, const Date&, std::vector< std::vector<MeteoData> >&,
                            const size_t&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GrassIO::writeMeteoData(const std::vector< std::vector<MeteoData> >&,
                             const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GrassIO::readPOI(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void GrassIO::write2DGrid(const Grid2DObject& grid_in, const std::string& name)
{
	fout.open(name.c_str());
	if (fout.fail()) {
		throw FileAccessException(name, AT);
	}

	Coords llcorner=grid_in.llcorner;
	//we want to make sure that we are using the provided projection parameters
	//so that we output is done in the same system as the inputs
	llcorner.setProj(coordout, coordoutparam);

	fout << setprecision(6) << fixed;

	try {
		fout << "north:" << (llcorner.getNorthing()+grid_in.cellsize*(double)grid_in.nrows) << "\n";
		fout << "south:" << llcorner.getNorthing() << "\n";
		fout << "east:"  << (llcorner.getEasting()+grid_in.cellsize*(double)grid_in.ncols)  << "\n";
		fout << "west:"  << llcorner.getEasting() << "\n";
		fout << "rows:"  << grid_in.nrows << "\n";
		fout << "cols:"  << grid_in.ncols << "\n";

		if(grid_in.nrows>0) {
			for (size_t kk=grid_in.nrows-1; kk < grid_in.nrows; kk--) {
				size_t ll = 0;
				for (ll=0; ll < (grid_in.ncols-1); ll++){
					if (grid_in(ll,kk) == IOUtils::nodata) {
						fout << "* ";
					} else {
						fout << grid_in(ll, kk) << " ";
					}
				}

				//The last value in a line does not have a trailing " "
				if (grid_in(ll,kk) == IOUtils::nodata) {
					fout << "*";
				} else {
					fout << grid_in(ll, kk);
				}
				fout << "\n";
			}
		}
	} catch(const std::exception& e) {
		cerr << "[E] " << AT << ": " << e.what() << std::endl;
		cleanup();
		throw;
	}

	cleanup();
}

void GrassIO::write2DGrid(const Grid2DObject&, const MeteoGrids::Parameters&, const Date&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

} //namespace
