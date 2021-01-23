cmake_minimum_required(VERSION 3.12)

option (WITH_CGAL "Include CGAL. Valgrind will puke with CGAL (rounding) so turn it off to use it." ON)
option (ARM64 "Compile for 64 bit ARM" OFF)

if (ARM64) 
	set (CMAKE_SYSTEM_NAME Linux)
	set (CMAKE_SYSTEM_PROCESSOR arm)
	set (CMAKE_C_COMPILER aarch64-linux-gnu-gcc-8)
	set (CMAKE_CXX_COMPILER aarch64-linux-gnu-g++-8)
endif (ARM64)

project (libgeo)

enable_language(C)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set (CMAKE_INSTALL_RPATH "/usr/local/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(${WIN32})
	message (STATUS "Building for Windows")
	set (CMAKE_CXX_FLAGS "/Zc:alignedNew /Zc:__cplusplus /std:c++17")	
	set (CMAKE_CXX_FLAGS_DEBUG "/Od /Zi /GS- /EHsc /guard:cf-")	
	add_definitions (-DBOOST_ALL_NO_LIB)
	set (BOOST_ROOT "C:\\dev\\boost_1_67_0")
	set (DEV_DIR "C:\\dev")
	set (GIT_DIR "${DEV_DIR}\\git")
	set (GIS_DIR "${DEV_DIR}\\OSGeo4W64")
	set (GEOS_LIBRARYDIR "${GIS_DIR}\\lib;${GIS_DIR}\\bin")
	set (GEOS_INCLUDEDIR "${GIS_DIR}\\include")
	set (EIGEN3_INCLUDE_DIR "${GIS_DIR}\\include\\eigen3")
	set (GDAL_INCLUDE_DIR "${GIS_DIR}\\include")
	set (GDAL_LIBRARYDIR "${GIS_DIR}\\lib;${GIS_DIR}\\bin")
	set (GIS_INCLUDE_DIR "${GIS_DIR}\\include")
	set (GIS_LIBRARYDIR "${GIS_DIR}\\lib;${GIS_DIR}\\bin")
	set (CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${GIS_DIR}\\bin;${DEV_DIR}\\Qt\\5.14.2\\msvc2017_64\\lib\\cmake")
	set (SQLITE_LIBRARY sqlite3_i)
	set (SPATIALITE_LIBRARY spatialite_i)
	set (GDAL_LIBRARY gdal_i)
	set (CGAL_INCLUDE_DIRS "C:\\dev\\CGAL-5.0.2\\include")
	set (GMP_LIBDIR "C:\\dev\\CGAL-5.0.2\\auxiliary\\gmp\\lib")
	set (GMP_LIBRARY libgmp-10)
	find_package (Boost 1.67)
else()
	message (STATUS "Building for Linux")
	set (CMAKE_CXX_COMPILER_LAUNCHER ccache)
	set (CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unknown-pragmas -fopenmp -frounding-math -D_GLIBCXX_PARALLEL" CACHE STRING "compile flags" FORCE)
	set (CMAKE_CXX_FLAGS_DEBUG "-g -O0" CACHE STRING "compile flags" FORCE)
	set (CMAKE_CXX_FLAGS_RELEASE "-O3" CACHE STRING "compile flags" FORCE)
	set (FS_LIB stdc++fs) # For filesystem features in c++17.
	find_package (Eigen3 REQUIRED)
	find_package (PDAL REQUIRED)
	find_package (OpenSSL REQUIRED)
	find_package (GDAL 2 REQUIRED)
	if (WITH_CGAL)
		find_package (CGAL REQUIRED)
		set (CGAL_LIBS "CGAL")
	else()
		message ("CGAL disabled.")
	endif (WITH_CGAL)
endif()

set (THREADS_PREFER_PTHREAD_FLAG ON)

find_package (Threads REQUIRED)

find_library (JSONCPP_LIBRARY jsoncpp)
#find_library (CRYPTO_LIBRARY crypto libeay32 PATHS ${GIS_LIBRARYDIR})
find_library (GEOS_LIBRARY geos_c PATHS ${GEOS_LIBRARYDIR})
find_library (PROJ_LIBRARY proj)

# Configure directories###################################################################################

include_directories(
	include
	ann/include 
	${OPENSSL_INCLUDE_DIR}  
	${GEOS_INCLUDEDIR} 
	${GIS_INCLUDEDIR} 
	${EIGEN3_INCLUDE_DIR} 
	${JSON_INCLUDEDIR}
	${GDAL_INCLUDE_DIR}
	${CGAL_INCLUDE_DIRS}
	${PDAL_INCLUDE_DIRS}
	${Boost_INCLUDE_DIRS}
)
	 
link_directories(./build/lib ./build/bin 
	${GDAL_LIBRARYDIR} 
	${GIS_LIBRARYDIR} 
	${PDAL_LIBRARY_DIRS}
	${GEOS_INCLUDEDIR}
	${PDAL_LIBRARY_DIRS}
	${GMP_LIBDIR}
	${Boost_LIBRARY_DIRS}
	${CGAL_DIR}
)

# Build support libraries ################################################################################

if (${WITH_FITPACK})
	message (STATUS "Building FitPack...")
	enable_language (Fortran)
	file (GLOB FSRC src/fitpack/*.f)
else ()
	file (GLOB FSRC src/geo.cpp src/fitpack/fitpack.cpp)
endif ()
add_library (fitpack_mod SHARED ${FSRC})
target_compile_features(fitpack_mod PUBLIC cxx_std_17)

set (geoutilsrc src/geo.cpp src/util.cpp src/md5.cpp)
add_library (geoutilobj OBJECT ${geoutilsrc})
target_compile_features(geoutilobj PUBLIC cxx_std_17)
set_target_properties (geoutilobj PROPERTIES POSITION_INDEPENDENT_CODE ON)
add_library(geoutil SHARED $<TARGET_OBJECTS:geoutilobj>)
target_link_libraries (geoutil ${FS_LIB} fitpack_mod ${GDAL_LIBRARY} ${CRYPTO_LIBRARY} ${SQLITE_LIBRARY} ${SPATIALITE_LIBRARY})

set (geodbsrc src/geo.cpp src/db.cpp)
add_library (geodbobj OBJECT ${geodbsrc})
target_compile_features(geodbobj PUBLIC cxx_std_17)
add_library (geodb SHARED $<TARGET_OBJECTS:geodbobj>)
set_target_properties (geodb PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries (geodb geoutil ${GDAL_LIBRARY})

set (geogridsrc src/geo.cpp src/grid.cpp)
add_library (geogridobj OBJECT ${geogridsrc})
target_compile_features(geogridobj PUBLIC cxx_std_17)
add_library (geogrid SHARED $<TARGET_OBJECTS:geogridobj>)
set_target_properties (geogrid PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries (geogrid geoutil geodb Threads::Threads ${GEOS_LIBRARY} ${GDAL_LIBRARY})

if (WITH_ANN)
	message (STATUS "Building ANN...")
	set (geoannsrc src/geo.cpp ann/src/ANN.cpp ann/src/brute.cpp 
		ann/src/kd_tree.cpp ann/src/kd_util.cpp 
		ann/src/kd_split.cpp ann/src/kd_dump.cpp ann/src/kd_search.cpp ann/src/kd_pr_search.cpp 
		ann/src/kd_fix_rad_search.cpp ann/src/bd_tree.cpp ann/src/bd_search.cpp ann/src/bd_pr_search.cpp 
		ann/src/bd_fix_rad_search.cpp ann/src/perf.cpp)
	add_library (geoannobj OBJECT ${geoannsrc})
	target_compile_features(geoannobj PUBLIC cxx_std_17)
	set_target_properties (geoannobj PROPERTIES POSITION_INDEPENDENT_CODE ON)
	target_compile_definitions (geoannobj PRIVATE DLL_EXPORTS=1)
	add_library (geoann SHARED $<TARGET_OBJECTS:geoannobj>)
	target_link_libraries (geoann geoutil)
	set (ANN_LIB geoann)
endif (WITH_ANN)

if (WITH_PC)
	message (STATUS "Building PC...")
	file(GLOB PC_COMPUTERS src/pc/computers/*.cpp)
	set (geopcsrc src/geo.cpp src/pointcloud.cpp src/pc/normalizer.cpp src/pc/rasterizer.cpp src/pc/filter.cpp ${PC_COMPUTERS})
	add_library (geopcobj OBJECT ${geopcsrc})
	target_compile_definitions (geopcobj PUBLIC DISABLE_CGAL)
	target_compile_features(geopcobj PUBLIC cxx_std_17)
	set_target_properties (geopcobj PROPERTIES POSITION_INDEPENDENT_CODE ON)
	target_compile_definitions (geopcobj PRIVATE DLL_EXPORTS=1)
	add_library (geopc SHARED $<TARGET_OBJECTS:geopcobj>)
	target_link_libraries (geopc geoutil geogrid ${CGAL_LIBS} ${ANN_LIB} ${PDAL_LIBRARIES} 
		${GMP_LIBRARY} ${Boost_LIBRARIES} ${GDAL_LIBRARIES})
	set (PC_LIB geopc)
endif (WITH_PC)

# Applications ##################################################################################################

if(${BUILD_APPS})
	message (STATUS "Building Apps...")
	
	# add_executable (pcperturb src/apps/pcperturb.cpp)
	# target_link_libraries (pcperturb geoutil geopc ${PDAL_LIBRARIES})

	# add_executable (pcregister src/apps/pcregister.cpp src/pc/trajectory.cpp)
	# target_link_libraries (pcregister geoutil geopc ${PDAL_LIBRARIES} ${PCL_LIBRARIES})

	# add_executable (geoinfo src/apps/geoinfo.cpp)
	# target_link_libraries (geoinfo geoutil ${GDAL_LIBRARY} ${PDAL_LIBRARIES} ${JSONCPP_LIBRARY})

	# add_executable (voidfill src/apps/voidfill.cpp)
	# target_link_libraries (voidfill geoutil geogrid ${GDAL_LIBRARY})

	# add_executable (datumshift src/apps/datumshift.cpp)
	# target_link_libraries (datumshift geoutil geogrid ${GDAL_LIBRARY} ${PROJ_LIBRARY})

	# add_executable (extmergesort src/apps/externalmergesort.cpp)
	# target_link_libraries (extmergesort geoutil geopc ${PDAL_LIBRARIES})

	# add_executable (groundassign src/apps/groundassign.cpp)
	# target_link_libraries (groundassign ${PDAL_LIBRARIES})

	add_executable (rastermatch src/apps/rastermatch.cpp)
	target_link_libraries (rastermatch geoutil ${ANN_LIB} geogrid ${GEOS_LIBRARY} ${GDAL_LIBRARY})

	# add_executable (raster_test src/raster.cpp)
	# target_link_libraries (raster_test geoutil ${ANN_LIB} geogrid ${GEOS_LIBRARY} ${GDAL_LIBRARY})

	# add_executable (cloudmatch src/apps/cloudmatch.cpp)
	# target_link_libraries (cloudmatch ${ANN_LIB} ${GDAL_LIBRARY} ${PDAL_LIBRARIES} fitpack_mod)

	# add_executable (rbf src/tests/rbf.cpp)
	# target_link_libraries (rbf)

	install(TARGETS ${ANN_LIB} ${PC_LIB} geogrid geoutil geodb fitpack_mod RUNTIME DESTINATION bin LIBRARY DESTINATION lib) 
	
else()

	install(TARGETS ${ANN_LIB} ${PC_LIB} geogrid geoutil geodb fitpack_mod RUNTIME DESTINATION bin LIBRARY DESTINATION lib) 

endif()


if(${BUILD_TESTS})

	#add_executable (libgeo_test src/tests/test.cpp)
	#target_link_libraries (libgeo_test geoutil geogrid geopc)

	#add_executable (bvspline_test src/tests/bvspline_test.cpp)
	#target_link_libraries (bvspline_test geoutil fitpack_mod)

endif()
