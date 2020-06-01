cmake_minimum_required(VERSION 3.12)

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
	set (DEV_DIR "C:\\dev")
	set (GIT_DIR "${DEV_DIR}\\git")
	set (GIS_DIR "${DEV_DIR}\\gisinternals")
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
	set (CGAL_DIR "C:\\dev\\CGAL-4.13")
	set (CMAKE_CXX_FLAGS "/Zc:alignedNew /Zc:__cplusplus /std:c++17")	
	set (CMAKE_CXX_FLAGS_DEBUG "/Od /O0 /Zi /GS- /guard:cf-")	
else()
	set (CMAKE_CXX_COMPILER_LAUNCHER ccache)
	set (CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unknown-pragmas -fopenmp -D_GLIBCXX_PARALLEL -std=c++17 -lstdc++fs")
	set (CMAKE_CXX_FLAGS_DEBUG "-g -O0")
	set (CMAKE_CXX_FLAGS_RELEASE "-O3")
	set (FS_LIB stdc++fs) # For filesystem features in c++17.
	find_package (Eigen3 REQUIRED)
endif()

set (THREADS_PREFER_PTHREAD_FLAG ON)

find_package (Threads REQUIRED)
find_package (GDAL 2 REQUIRED)

if (${WITH_CGAL})
	find_package (CGAL REQUIRED)
endif()

if (${WITH_PDAL})
	find_package (PDAL REQUIRED)
endif()

if (${WITH_OPENSSL})
	find_package (OpenSSL REQUIRED)
endif()

find_library (JSONCPP_LIBRARY jsoncpp)
#find_library (CRYPTO_LIBRARY crypto libeay32 PATHS ${GIS_LIBRARYDIR})
find_library (GEOS_LIBRARY geos_c PATHS ${GEOS_LIBRARYDIR})
find_library (PROJ_LIBRARY proj)

if (${WITH_LIBLAS})
	find_library (LIBLAS_LIBRARY las)
endif()

# Configure directories###################################################################################

include_directories(BEFORE SYSTEM 
	include
	ann/include 
	${OPENSSL_INCLUDE_DIR}  
	${GEOS_INCLUDEDIR} 
	${GIS_INCLUDEDIR} 
	${EIGEN3_INCLUDE_DIR} 
	${JSON_INCLUDEDIR}
	${GDAL_INCLUDE_DIR}
)
	 
link_directories(./build/lib ./build/bin 
	${GDAL_LIBRARYDIR} 
	${GIS_LIBRARYDIR} 
	${PDAL_LIBRARY_DIRS}
	${GEOS_INCLUDEDIR}
)

# Build support libraries ################################################################################

if (${WITH_FITPACK})
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

if(${WITH_ANN})
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
endif()

if (${WITH_PC})
	file(GLOB PC_COMPUTERS src/pc/computers/*.cpp)
	add_library (geopc SHARED src/geo.cpp src/pointcloud.cpp src/pc/normalizer.cpp src/pc/rasterizer.cpp src/pc/filter.cpp ${PC_COMPUTERS})
	target_compile_features(geopc PUBLIC cxx_std_17)
	target_link_libraries (geopc geoutil geogrid ${ANN_LIB} ${LIBLAS_LIBRARY} CGAL CGAL_Core gmp)
	set (PC_LIB geopc)
endif()

# Applications ##################################################################################################

if(${BUILD_APPS})

	# add_executable (pc2tile src/apps/pc2tile.cpp)
	# target_link_libraries (pc2tile geoutil geopc ${LIBLAS_LIBRARY} ${GEOS_LIBRARY})

	# add_executable (pcperturb src/apps/pcperturb.cpp)
	# target_link_libraries (pcperturb geoutil geopc ${LIBLAS_LIBRARY})

	# add_executable (pcregister src/apps/pcregister.cpp src/pc/trajectory.cpp)
	# target_link_libraries (pcregister geoutil geopc ${LIBLAS_LIBRARY} ${PCL_LIBRARIES})

	# add_executable (geoinfo src/apps/geoinfo.cpp)
	# target_link_libraries (geoinfo geoutil ${GDAL_LIBRARY} ${LIBLAS_LIBRARY} ${JSONCPP_LIBRARY})

	add_executable (voidfill src/apps/voidfill.cpp)
	target_link_libraries (voidfill geoutil geogrid ${GDAL_LIBRARY})

	# add_executable (datumshift src/apps/datumshift.cpp)
	# target_link_libraries (datumshift geoutil geogrid ${GDAL_LIBRARY} ${PROJ_LIBRARY})

	# add_executable (extmergesort src/apps/externalmergesort.cpp)
	# target_link_libraries (extmergesort geoutil geopc ${LIBLAS_LIBRARY})

	# add_executable (groundassign src/apps/groundassign.cpp)
	# target_link_libraries (groundassign ${LIBLAS_LIBRARY})

	# add_executable (rastermatch src/apps/rastermatch.cpp)
	# target_link_libraries (rastermatch geoutil ${ANN_LIB} geogrid ${GEOS_LIBRARY} ${GDAL_LIBRARY})

	# add_executable (raster_test src/raster.cpp)
	# target_link_libraries (raster_test geoutil ${ANN_LIB} geogrid ${GEOS_LIBRARY} ${GDAL_LIBRARY})

	# add_executable (cloudmatch src/apps/cloudmatch.cpp)
	# target_link_libraries (cloudmatch ${ANN_LIB} ${GDAL_LIBRARY} ${PDAL_LIBRARIES} fitpack_mod)

	# add_executable (rbf src/tests/rbf.cpp)
	# target_link_libraries (rbf)

	install(TARGETS voidfill ${ANN_LIB} ${PC_LIB} geogrid geoutil geodb fitpack_mod RUNTIME DESTINATION bin LIBRARY DESTINATION lib) 
	
else()

	install(TARGETS ${ANN_LIB} ${PC_LIB} geogrid geoutil geodb fitpack_mod RUNTIME DESTINATION bin LIBRARY DESTINATION lib) 

endif()


if(${BUILD_TESTS})

	#add_executable (libgeo_test src/tests/test.cpp)
	#target_link_libraries (libgeo_test geoutil geogrid geopc)

	#add_executable (bvspline_test src/tests/bvspline_test.cpp)
	#target_link_libraries (bvspline_test geoutil fitpack_mod)

endif()
