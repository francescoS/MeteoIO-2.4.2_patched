/***********************************************************************************/
/*  Copyright 2014 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef __LIBNCPP_H__
#define __LIBNCPP_H__

#include <meteoio/IOExceptions.h>
#include <meteoio/IOUtils.h>

#include <netcdf.h>
#include <string>
#include <vector>

namespace ncpp {

	//Opening, creating, closing dataset
	void open_file(const std::string& filename, const int& omode, int& ncid);
	void create_file(const std::string& filename, const int& cmode, int& ncid);
	void start_definitions(const std::string& filename, const int& ncid);
	void end_definitions(const std::string& filename, const int& ncid);
	void close_file(const std::string& filename, const int& ncid);

	//Adding variables
	void add_0D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, int& varid);
	void add_1D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid, int& varid);
	void add_2D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid1, const int& dimid2, int& varid);
	void add_3D_variable(const int& ncid, const std::string& varname, const nc_type& xtype, const int& dimid_record,
	                     const int& dimid1, const int& dimid2, int& varid);

	//Adding attributes
	void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const std::string& attr_value);
	void add_attribute(const int& ncid, const int& varid, const std::string& attr_name, const double& attr_value);
	void get_attribute(const int& ncid, const std::string& varname, const int& varid, const std::string& attr_name, std::string& attr_value);
	void get_attribute(const int& ncid, const std::string& varname, const int& varid, const std::string& attr_name, double& attr_value);
	bool check_attribute(const int& ncid, const int& varid, const std::string& attr_name);

	//Adding dimensions
	void add_dimension(const int& ncid, const std::string& dimname, const size_t& length, int& dimid);


	//Reading data from NetCDF file
	void read_data(const int& ncid, const std::string& varname, const int& varid,
	               const size_t& pos, const size_t& latlen, const size_t& lonlen, double*& data);
	void read_data_2D(const int& ncid, const std::string& varname, const int& varid,
	                  const size_t& record, const size_t& count, const size_t& length, double*& data);
	void read_value(const int& ncid, const std::string& varname, const int& varid, double& data);
	void read_value(const int& ncid, const std::string& varname, const int& varid, const size_t& pos, double& data);
	void read_data(const int& ncid, const std::string& varname, const int& varid, double*& data);

	//Writing data to NetCDF file
	void write_data(const int& ncid, const std::string& varname, const int& varid, const double * const data);
	void write_data(const int& ncid, const std::string& varname, const int& varid, const size_t& nrows, const size_t& ncols,
	                const size_t& pos_start, const double * const data);

	//Dealing with variables that have dimension NC_UNLIMITED
	size_t find_record(const int& ncid, const std::string& varname, const int& varid, const double& data);
	size_t add_record(const int& ncid, const std::string& varname, const int& varid, const double& data);
	void write_record(const int& ncid, const std::string& varname, const int& varid, const size_t& pos,
	                  const size_t& length, const double * const data);

	//Dealing with variables and dimensions
	bool check_dim_var(const int& ncid, const std::string& dimname);
	bool check_variable(const int& ncid, const std::string& varname);
	void get_variable(const int& ncid, const std::string& varname, int& varid);
	void get_variables(const int& ncid, const std::vector<std::string>& dimensions, std::vector<std::string>& variables);
	bool check_dimensions(const int& ncid, const std::string& varname, const int& varid, const std::vector<std::string>& names);
	void get_dimension(const int& ncid, const std::string& dimname, int& dimid);
	void get_dimension(const int& ncid, const std::string& dimname, int& dimid, size_t& dimlen);
	void get_dimension(const int& ncid, const std::string& varname, const int& varid,
	                   std::vector<int>& dimid, std::vector<int>& dim_varid, std::vector<std::string>& dimname, std::vector<size_t>& dimlen);

} // end namespace

#endif
