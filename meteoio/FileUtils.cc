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

#include <algorithm>
#include <iostream>
#include <cstdio>
#include <fstream>

#if defined _WIN32 || defined __MINGW32__
	#include <windows.h>
	#include "Shlwapi.h"
	//removing two macros defined in windows.h
	#undef max
	#undef min
#else
	#include <dirent.h>
	#include <sys/stat.h>
#endif

#include <meteoio/FileUtils.h>
#include <meteoio/IOUtils.h>

namespace mio {
namespace IOUtils {

void copy_file(const std::string& src, const std::string& dest)
{
	if (src == dest) return; //copying to the same file doesn't make sense, but is no crime either

	std::ifstream fin(src.c_str(), std::ios::binary);
	if (fin.fail()) throw FileAccessException(src, AT);

	std::ofstream fout(dest.c_str(), std::ios::binary);
	if (fout.fail()) {
		fin.close();
		throw FileAccessException(dest, AT);
	}

	fout << fin.rdbuf();

	fin.close();
	fout.close();
}

std::string cleanPath(std::string in_path, const bool& resolve)
{
	if(!resolve) { //do not resolve links, relative paths, etc
		std::replace(in_path.begin(), in_path.end(), '\\', '/');
		return in_path;
	} else {
	#if defined _WIN32 || defined __MINGW32__
		//if this would not suffice, see http://pdh11.blogspot.ch/2009/05/pathcanonicalize-versus-what-it-says-on.html
		char **ptr = NULL;
		char *out_buff = (char*)calloc(MAX_PATH, sizeof(char));
		const DWORD status = GetFullPathName(in_path.c_str(), MAX_PATH, out_buff, ptr);
		if(status!=0 && status<=MAX_PATH) in_path = out_buff;
		free(out_buff);

		std::replace(in_path.begin(), in_path.end(), '\\', '/');
		return in_path;
	#else //POSIX
		std::replace(in_path.begin(), in_path.end(), '\\', '/');

		char *real_path = realpath(in_path.c_str(), NULL); //POSIX
		if(real_path!=NULL) {
			const std::string tmp(real_path);
			free(real_path);
			return tmp;
		} else
			return in_path; //something failed in realpath, keep it as it is
	#endif
	}
}

std::string getExtension(const std::string& filename)
{
	const size_t start_basename = filename.find_last_of("/\\"); //we will skip the path
	const size_t startpos = filename.find_last_of('.');
	if( startpos==std::string::npos ) return std::string();
	if( start_basename!=std::string::npos && startpos<start_basename ) return std::string();

	const std::string whitespaces(" \t\f\v\n\r");
	const size_t endpos = filename.find_last_not_of(whitespaces); // Find the first character position from reverse af

	return filename.substr(startpos+1, endpos-startpos);
}

std::string removeExtension(const std::string& filename)
{
	const size_t start_basename = filename.find_last_of("/\\"); //we will skip the path
	const size_t startpos = filename.find_last_of('.');
	if( startpos==std::string::npos ) return filename;
	if( start_basename!=std::string::npos && startpos<start_basename ) return filename;

	return filename.substr(0, startpos);
}

std::string getPath(const std::string& filename, const bool& resolve)
{
	const std::string clean_filename = cleanPath(filename, resolve);
	const size_t end_path = clean_filename.find_last_of("/");
	if(end_path!=std::string::npos) {
		return clean_filename.substr(0, end_path);
	} else {
		return cleanPath("./", resolve);
	}
}

std::string getFilename(const std::string& path)
{
	const size_t start_basename = path.find_last_of("/\\");
	if(start_basename!=std::string::npos)
		return path.substr(start_basename+1, std::string::npos);
	else
		return path;
}

bool validFileName(const std::string& filename)
{
	const size_t startpos = filename.find_first_not_of(" \t\n"); // Find the first character position after excluding leading blank spaces
	if((startpos!=0) || (filename==".") || (filename=="..")) {
		return false;
	}

	return true;
}


bool isAbsolutePath(const std::string& in_path)
{
#if defined _WIN32 || defined __MINGW32__ || defined __CYGWIN__
	return (in_path.size()>=2 && in_path[1]==':');
#else
	return (in_path.size()>=1 && in_path[0]=='/');
#endif
}

#if defined _WIN32 || defined __MINGW32__
bool fileExists(const std::string& filename)
{
	return ( GetFileAttributes( filename.c_str() ) != INVALID_FILE_ATTRIBUTES );
}

void readDirectory(const std::string& path, std::list<std::string>& dirlist, const std::string& pattern)
{
	const size_t path_length = path.length();
	if (path_length > (MAX_PATH - 3)) {
		std::cerr << "Path " << path << "is too long (" << path_length << " characters)" << std::endl;
		throw FileAccessException("Error opening directory " + path, AT);
	}

	const std::string filepath = path+"\\"+pattern;
	WIN32_FIND_DATA ffd;
	const HANDLE hFind = FindFirstFile(filepath.c_str(), &ffd);
	if (INVALID_HANDLE_VALUE == hFind) {
		throw FileAccessException("Error opening directory " + path, AT);
	}

	do {
		if (ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
			//this is a directory -> do nothing
		} else {
			const std::string filename(ffd.cFileName);
			dirlist.push_back(filename);
		}
	}
	while (FindNextFile(hFind, &ffd) != 0);

	const DWORD dwError = GetLastError();
	if (dwError != ERROR_NO_MORE_FILES) {
		throw FileAccessException("Error listing files in directory " + path, AT);
	}

	FindClose(hFind);
}
#else
bool fileExists(const std::string& filename)
{
	struct stat buffer ;

	if ((stat( filename.c_str(), &buffer))==0) {//File exists if stat returns 0
		return true ;
	}

	return false;
}

void readDirectory(const std::string& path, std::list<std::string>& dirlist, const std::string& pattern)
{
	DIR *dp = opendir(path.c_str());
	if(dp == NULL) {
		throw FileAccessException("Error opening directory " + path, AT);
	}

	struct dirent *dirp;
	while ((dirp = readdir(dp)) != NULL) {
		const std::string tmp(dirp->d_name);
		if( tmp.compare(".")!=0 && tmp.compare("..")!=0 ) { //skip "." and ".."
			if (pattern.empty()) {
				dirlist.push_back(tmp);
			} else {
				const size_t pos = tmp.find(pattern);
				if (pos!=std::string::npos) {
					dirlist.push_back(tmp);
				}
			}
		}
	}
	closedir(dp);
}
#endif

char getEoln(std::istream& fin)
{
	std::streambuf* pbuf;
	char tmp = '0';
	int chars = 0;

	const std::streampos file_start = fin.tellg();

	do {
		fin.get(tmp);
		chars++;

		if ((tmp == '\r') || (tmp == '\n')) {
			char peekc = tmp;
			while ((!fin.eof() && ((peekc=='\r') || (peekc=='\n')))) {
				tmp = peekc;
				fin.get(peekc);
				chars++;
			}
			pbuf = fin.rdbuf();
			pbuf->pubseekpos(file_start); //rewind
			fin.clear(); //reset eof flag, etc
			return tmp;
		}
	} while ((chars < 3000) && (!fin.eof()));

	pbuf = fin.rdbuf();
	pbuf->pubseekpos(file_start); //rewind
	fin.clear(); //reset eof flag, etc

	return '\n';
}

void skipLines(std::istream& fin, const size_t& nbLines, const char& eoln)
{
	std::string dummy;
	for (size_t ii=0; ii<nbLines; ii++) {
		if(!getline(fin, dummy, eoln)) {
			throw InvalidFormatException("Premature EOF while skipping lines", AT);
		}
	}
}


//below, the file indexer implementation
void FileIndexer::setIndex(const Date& i_date, const std::streampos& i_pos)
{
	const file_index elem(i_date, i_pos);

	//check if we can simply append the new index
	if(vecIndex.empty() || elem>vecIndex.back()) {
		vecIndex.push_back(elem);
		return;
	}

	//look for the proper position for insertion of the new index
	const std::vector< struct file_index >::iterator it = std::upper_bound(vecIndex.begin(), vecIndex.end(), elem);
	if(it>vecIndex.begin() && (it-1)->date!=elem.date) { //check that we don't try to insert a duplicate
		vecIndex.insert(it, elem); //insertion is at the proper place -> remains ordered
		return;
	}
}

void FileIndexer::setIndex(const std::string& i_date, const std::streampos& i_pos)
{
	Date tmpdate;
	convertString(tmpdate, i_date, 0.);
	setIndex(tmpdate, i_pos);
}

void FileIndexer::setIndex(const double& i_date, const std::streampos& i_pos)
{
	const Date tmpdate(i_date, 0.);
	setIndex(tmpdate, i_pos);
}

std::streampos FileIndexer::getIndex(const Date& i_date) const
{
	const size_t foundIdx = binarySearch(i_date);
	if(foundIdx==static_cast<size_t>(-1)) return static_cast<std::streampos>(-1);
	else return vecIndex[foundIdx].pos;
}

std::streampos FileIndexer::getIndex(const std::string& i_date) const
{
	Date tmpdate;
	convertString(tmpdate, i_date, 0.);
	return getIndex(tmpdate);
}

std::streampos FileIndexer::getIndex(const double& i_date) const
{
	const Date tmpdate(i_date, 0.);
	return getIndex(tmpdate);
}

size_t FileIndexer::binarySearch(const Date& soughtdate) const
{//perform binary search, return the first element that is GREATER than the provided value
	if(vecIndex.empty()) return static_cast<size_t>(-1);
	if(soughtdate<vecIndex.front().date) return static_cast<size_t>(-1);
	if(soughtdate>=vecIndex.back().date) return vecIndex.size()-1;

	const file_index elem(soughtdate, 0);
	//returns the first element that is GREATER than the provided value
	const std::vector< struct file_index >::const_iterator it = std::upper_bound(vecIndex.begin(), vecIndex.end(), elem);
	if(it>vecIndex.begin()) return it-vecIndex.begin()-1;
	else return static_cast<size_t>(-1);
}

const std::string FileIndexer::toString() const
{
	std::ostringstream os;
	os << "<FileIndexer>\n";
	for(size_t ii=0; ii<vecIndex.size(); ii++)
		os << "\t" << "[" << ii << "] - " << vecIndex[ii].date.toString(Date::ISO) << " -> #" << std::hex << vecIndex[ii].pos << std::dec << "\n";
	os << "</FileIndexer>\n";
	return os.str();
}

} //end namespace IOUtils
} //end namespace mio
