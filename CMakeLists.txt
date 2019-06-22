cmake_minimum_required(VERSION 3.1)

project (libgeo)

enable_language(C)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -D_GLIBCXX_PARALLEL")
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set (CMAKE_INSTALL_RPATH "/usr/local/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set (Boost_DEBUG OFF)
set (Boost_USE_STATIC_LIBS OFF)
set (Boost_USE_STATIC_RUNTIME OFF)
set (Boost_USE_MULTITHREADED ON)
set (Boost_NO_SYSTEM_PATHS OFF)

if (${APPLE})
	set (CMAKE_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
elseif(${WIN32})
	set (DEV_DIR "C:\\dev")
	set (Boost_COMPILER -vc141)
	set (BOOST_ROOT "${DEV_DIR}\\boost_1_67_0")
	set (GIT_DIR "${DEV_DIR}\\git")
	set (GIS_DIR "${DEV_DIR}\\gisinternals")
	set (GEOS_LIBRARYDIR "${GIS_DIR}\\lib;${GIS_DIR}\\bin") #set (GEOS_LIBRARYDIR "${DEV_DIR}\\geos\\lib;${DEV_DIR}\\geos\\bin")
	set (GEOS_INCLUDEDIR "${GIS_DIR}\\include") #set (GEOS_INCLUDEDIR "${DEV_DIR}\\geos\\include")
	set (EIGEN3_INCLUDE_DIR "${DEV_DIR}\\eigen3")
	set (GDAL_INCLUDE_DIR "${GIS_DIR}\\include")
	set (GDAL_LIBRARYDIR "${GIS_DIR}\\lib;${GIS_DIR}\\bin")
	set (GIS_INCLUDE_DIR "${GIS_DIR}\\include")
	set (GIS_LIBRARYDIR "${GIS_DIR}\\lib;${GIS_DIR}\\bin")
	set (CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${GIS_DIR}\\bin;${DEV_DIR}\\Qt\\5.11.0\\msvc2017_64\\lib\\cmake")
	set (SQLITE_LIBRARY sqlite3_i)
	set (SPATIALITE_LIBRARY spatialite_i)
	set (GDAL_LIBRARY gdal_i)
	set (CGAL_DIR "C:\\dev\\CGAL-4.13")
else()
	set (CMAKE_CXX_COMPILER_LAUNCHER ccache)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas -fopenmp")
endif()

set (THREADS_PREFER_PTHREAD_FLAG ON)

find_package (Threads REQUIRED)
find_package (GDAL 2 REQUIRED)
find_package (CGAL REQUIRED)
find_package (PDAL REQUIRED)
find_package (OpenSSL REQUIRED)
find_package (Eigen3 REQUIRED)
find_package (Boost 1.58 COMPONENTS system filesystem date_time REQUIRED)

find_library (LIBLAS_LIBRARY las)
find_library (JSONCPP_LIBRARY jsoncpp)
find_library (CRYPTO_LIBRARY crypto libeay32 PATHS ${GIS_LIBRARYDIR})
find_library (GEOS_LIBRARY geos_c PATHS ${GEOS_LIBRARYDIR})
find_library (PROJ_LIBRARY proj)


# Configure directories###################################################################################

include_directories(BEFORE SYSTEM ./include ./ann/include 
	${OPENSSL_INCLUDE_DIR}  
	${GEOS_INCLUDEDIR} ${GIS_INCLUDEDIR} ${EIGEN3_INCLUDE_DIR} ${JSON_INCLUDEDIR}
	${GDAL_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
	 
link_directories(./build/lib ./build/bin ${GDAL_LIBRARYDIR} 
	${GIS_LIBRARYDIR} ${Boost_LIBRARY_DIRS} ${PDAL_LIBRARY_DIRS})

# Build support libraries ################################################################################

set (geoutilsrc src/geo.cpp src/util.cpp src/md5.cpp)
add_library (geoutilobj OBJECT ${geoutilsrc})
set_target_properties (geoutilobj PROPERTIES POSITION_INDEPENDENT_CODE ON)
add_library(geoutil SHARED $<TARGET_OBJECTS:geoutilobj>)
target_link_libraries (geoutil ${Boost_LIBRARIES} ${GDAL_LIBRARY} ${CRYPTO_LIBRARY} ${SQLITE_LIBRARY} ${SPATIALITE_LIBRARY})
# is rt needed here doesn't work on windows.

set (geodbsrc src/geo.cpp src/db.cpp)
add_library (geodbobj OBJECT ${geodbsrc})
set_target_properties (geodbobj PROPERTIES POSITION_INDEPENDENT_CODE ON)
add_library (geodb SHARED $<TARGET_OBJECTS:geodbobj>)
target_link_libraries (geodb geoutil ${GDAL_LIBRARY})

set (geogridsrc src/geo.cpp src/grid.cpp)
add_library (geogridobj OBJECT ${geogridsrc})
set_target_properties (geogridobj PROPERTIES POSITION_INDEPENDENT_CODE ON)
add_library (geogrid SHARED $<TARGET_OBJECTS:geogridobj>)
target_link_libraries (geogrid geoutil geodb Threads::Threads ${GEOS_LIBRARY} ${GDAL_LIBRARY})

file(GLOB PC_COMPUTERS src/pc/computers/*.cpp)
add_library (geopc SHARED src/geo.cpp src/pointcloud.cpp src/pc/normalizer.cpp src/pc/rasterizer.cpp src/pc/filter.cpp ${PC_COMPUTERS})
target_link_libraries (geopc geoutil geogrid geoann ${LIBLAS_LIBRARY} CGAL CGAL_Core gmp)

set (geoannsrc src/geo.cpp ann/src/ANN.cpp ann/src/brute.cpp 
	ann/src/kd_tree.cpp ann/src/kd_util.cpp 
	ann/src/kd_split.cpp ann/src/kd_dump.cpp ann/src/kd_search.cpp ann/src/kd_pr_search.cpp 
	ann/src/kd_fix_rad_search.cpp ann/src/bd_tree.cpp ann/src/bd_search.cpp ann/src/bd_pr_search.cpp 
	ann/src/bd_fix_rad_search.cpp ann/src/perf.cpp)
add_library (geoannobj OBJECT ${geoannsrc})
set_target_properties (geoannobj PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_definitions (geoannobj PRIVATE DLL_EXPORTS=1)
add_library (geoann SHARED $<TARGET_OBJECTS:geoannobj>)
target_link_libraries (geoann geoutil)

# Applications ##################################################################################################

if(${BUILD_APPS})

	#enable_language (Fortran)
	#file (GLOB FSRC src/fitpack/*.f)
	#add_library (fitpack_mod SHARED ${FSRC}) 

	# add_executable (polygonize src/apps/polygonize.cpp)
	# target_link_libraries (polygonize geoutil geogrid ${Boost_LIBRARIES})

	add_executable (pc2grid src/apps/pc2grid.cpp)
	target_link_libraries (pc2grid geoann geoutil geogrid geopc geoann ${GEOS_LIBRARY})

	# add_executable (pc2tile src/apps/pc2tile.cpp)
	# target_link_libraries (pc2tile geoutil geopc ${Boost_LIBRARIES} ${LIBLAS_LIBRARY} ${GEOS_LIBRARY})

	add_executable (pcnorm src/apps/pcnorm.cpp)
	target_link_libraries (pcnorm ${LIBLAS_LIBRARY})

	# add_executable (pcperturb src/apps/pcperturb.cpp)
	# target_link_libraries (pcperturb geoutil geopc ${Boost_LIBRARIES} ${LIBLAS_LIBRARY})

	# add_executable (pcregister src/apps/pcregister.cpp src/pc/trajectory.cpp)
	# target_link_libraries (pcregister geoutil geopc ${Boost_LIBRARIES} ${LIBLAS_LIBRARY} ${PCL_LIBRARIES})

	# add_executable (geoinfo src/apps/geoinfo.cpp)
	# target_link_libraries (geoinfo geoutil ${GDAL_LIBRARY} ${LIBLAS_LIBRARY} ${JSONCPP_LIBRARY})

	#add_executable (voidfill src/apps/voidfill.cpp)
	#target_link_libraries (voidfill geoutil geogrid ${GDAL_LIBRARY})

	# add_executable (datumshift src/apps/datumshift.cpp)
	# target_link_libraries (datumshift geoutil geogrid ${GDAL_LIBRARY} ${PROJ_LIBRARY})

	# add_executable (extmergesort src/apps/externalmergesort.cpp)
	# target_link_libraries (extmergesort geoutil geopc ${LIBLAS_LIBRARY})

	# add_executable (sbet2csv src/apps/sbet2csv.cpp)
	# target_link_libraries (sbet2csv geoutil ${PROJ_LIBRARY})

	# add_executable (groundassign src/apps/groundassign.cpp)
	# target_link_libraries (groundassign ${LIBLAS_LIBRARY})

	#add_executable (rastermatch src/apps/rastermatch.cpp)
	#target_link_libraries (rastermatch geoutil geoann geogrid ${GEOS_LIBRARY} ${GDAL_LIBRARY})

	#add_executable (raster_test src/raster.cpp)
	#target_link_libraries (raster_test geoutil geoann geogrid ${GEOS_LIBRARY} ${GDAL_LIBRARY})

	#add_executable (cloudmatch src/apps/cloudmatch.cpp)
	#target_link_libraries (cloudmatch geoann ${Boost_LIBRARIES} ${GDAL_LIBRARY} ${PDAL_LIBRARIES} fitpack_mod)

	# add_executable (rbf src/tests/rbf.cpp)
	# target_link_libraries (rbf)

	# install(TARGETS pc2grid pc2tile rastermatch geopc geogrid geoinfo geoutil geoann RUNTIME DESTINATION bin LIBRARY DESTINATION lib) 
	
	install(TARGETS pc2grid pcnorm geoann geopc geogrid geoutil geodb RUNTIME DESTINATION bin LIBRARY DESTINATION lib) 
	
endif()

