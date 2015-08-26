CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../../tools/cmake/")

PROJECT(meteoio_examples)

FIND_PACKAGE(MeteoIO)
INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${METEOIO_INCLUDE_DIR})

OPTION(BUILD_SHARED_LIBS "Set to ON to build dynamic libraries" ON)

IF(NOT BUILD_SHARED_LIBS)
            FIND_PATH(LIBPROJ_INCLUDE_DIR "proj_api.h")
            FIND_LIBRARY(LIBPROJ proj)
            INCLUDE_DIRECTORIES(${LIBPROJ_INCLUDE_DIR})
            IF(MSVC)
                    ADD_DEFINITIONS(/DPROJ4) #it looks like some VC++ versions don't support -D syntax
            ELSE(MSVC)
                    ADD_DEFINITIONS(-DPROJ4)
            ENDIF(MSVC)
ENDIF(NOT BUILD_SHARED_LIBS)

#####################################################################
# add all the executables
ADD_EXECUTABLE(2D_interpolations 2D_interpolations.cc)
ADD_EXECUTABLE(coordinates coordinates.cc)
ADD_EXECUTABLE(data_converter data_converter.cc)
ADD_EXECUTABLE(dem_reading dem_reading.cc)
ADD_EXECUTABLE(grid2d_reading grid2d_reading.cc)
ADD_EXECUTABLE(matrix matrix.cc)
ADD_EXECUTABLE(meteo_reading meteo_reading.cc)
ADD_EXECUTABLE(sun sun.cc)
ADD_EXECUTABLE(time time.cc)

IF(BUILD_SHARED_LIBS)
    TARGET_LINK_LIBRARIES(2D_interpolations ${METEOIO_LIBRARIES})
    TARGET_LINK_LIBRARIES(coordinates ${METEOIO_LIBRARIES})
    TARGET_LINK_LIBRARIES(data_converter ${METEOIO_LIBRARIES})
    TARGET_LINK_LIBRARIES(dem_reading ${METEOIO_LIBRARIES})
    TARGET_LINK_LIBRARIES(grid2d_reading ${METEOIO_LIBRARIES})
    TARGET_LINK_LIBRARIES(matrix ${METEOIO_LIBRARIES})
    TARGET_LINK_LIBRARIES(meteo_reading ${METEOIO_LIBRARIES})
    TARGET_LINK_LIBRARIES(sun ${METEOIO_LIBRARIES})
    TARGET_LINK_LIBRARIES(time ${METEOIO_LIBRARIES})
ELSE(BUILD_SHARED_LIBS)
    TARGET_LINK_LIBRARIES(2D_interpolations ${METEOIO_LIBRARIES} ${LIBPROJ})
    TARGET_LINK_LIBRARIES(coordinates ${METEOIO_LIBRARIES} ${LIBPROJ})
    TARGET_LINK_LIBRARIES(data_converter ${METEOIO_LIBRARIES} ${LIBPROJ})
    TARGET_LINK_LIBRARIES(dem_reading ${METEOIO_LIBRARIES} ${LIBPROJ})
    TARGET_LINK_LIBRARIES(grid2d_reading ${METEOIO_LIBRARIES} ${LIBPROJ})
    TARGET_LINK_LIBRARIES(matrix ${METEOIO_LIBRARIES} ${LIBPROJ})
    TARGET_LINK_LIBRARIES(meteo_reading ${METEOIO_LIBRARIES} ${LIBPROJ})
    TARGET_LINK_LIBRARIES(sun ${METEOIO_LIBRARIES} ${LIBPROJ})
    TARGET_LINK_LIBRARIES(time ${METEOIO_LIBRARIES} ${LIBPROJ})
ENDIF(BUILD_SHARED_LIBS)

#####################################################################
# add a distclean target
ADD_CUSTOM_TARGET(distclean make clean
	COMMAND cmake -E remove *~ *.o
	COMMAND cmake -E remove CMakeCache.txt cmake_install.cmake
	COMMAND cmake -E remove_directory CMakeFiles
	)