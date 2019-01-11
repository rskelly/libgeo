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
set (CMAKE_INSTALL_RPATH "/opt/lib;/usr/local/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set (Boost_DEBUG OFF)
if(WIN32)
	set (Boost_USE_STATIC_LIBS ON)
else()
	set (Boost_USE_STATIC_LIBS OFF)
endif()
set (Boost_USE_STATIC_RUNTIME OFF)
set (Boost_USE_MULTITHREADED ON)
set (Boost_NO_SYSTEM_PATHS ON)

if (${APPLE})
	set (CMAKE_C_COMPILER /usr/local/bin/clang-omp CACHE STRING "C compiler" FORCE)
	set (CMAKE_CXX_COMPILER /usr/local/bin/clang-omp++ CACHE STRING "C++ compiler" FORCE)
	set (CMAKE_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
elseif(${WIN32})
	set (DEV_DIR "C:\\dev")
	set (Boost_COMPILER -vc141)
	set (BOOST_ROOT "${DEV_DIR}\\boost_1_67_0")
	set (GIT_DIR "${DEV_DIR}\\git")
	set (GIS_DIR "${DEV_DIR}\\gisinternals")
	set (GEOS_LIBRARYDIR "${DEV_DIR}\\geos\\lib;${DEV_DIR}\\geos\\bin")
	set (GEOS_INCLUDEDIR "${DEV_DIR}\\geos\\include")
	set (EIGEN_INCLUDEDIR "${DEV_DIR}\\eigen3")
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
	set (GEOS_INCLUDEDIR /opt/include)
	set (GEOS_LIBRARYDIR /opt/lib)
	set (EIGEN_INCLUDEDIR "/usr/include/eigen3")
	set (BOOST_INCLUDEDIR "/usr/include")
	set (BOOST_LIBRARYDIR "/usr/lib/x86_64-linux-gnu")
	set (JSON_INCLUDEDIR "/usr/include/jsoncpp")
endif()

set (THREADS_PREFER_PTHREAD_FLAG ON)

find_package (Threads REQUIRED)
find_package (Boost 1.58 COMPONENTS system filesystem date_time REQUIRED)
find_package (GDAL 2.2 REQUIRED)
find_package (CGAL REQUIRED)

find_library (LIBLAS_LIBRARY las)
find_library (JSONCPP_LIBRARY jsoncpp)
find_library (CRYPTO_LIBRARY crypto libeay32 PATHS ${GIS_LIBRARYDIR})
find_library (GEOS_LIBRARY geos PATHS ${GEOS_LIBRARYDIR})
find_library (PROJ_LIBRARY proj)

# Configure directories###################################################################################

include_directories(BEFORE SYSTEM ./include ./ann/include /opt/include
	${GEOS_INCLUDEDIR} ${GIS_INCLUDEDIR} ${EIGEN_INCLUDEDIR} ${JSON_INCLUDEDIR}
	${GDAL_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
	 
link_directories(./build/lib ./build/bin /opt/lib ${GDAL_LIBRARYDIR} 
	${GIS_LIBRARYDIR} ${Boost_LIBRARY_DIRS})

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

set (georastersrc src/geo.cpp src/raster.cpp src/rast/polygonize.cpp)
add_library (georasterobj OBJECT ${georastersrc})
set_target_properties (georasterobj PROPERTIES POSITION_INDEPENDENT_CODE ON)
add_library (georaster SHARED $<TARGET_OBJECTS:georasterobj>)
target_link_libraries (georaster geoutil geodb Threads::Threads ${GEOS_LIBRARY} ${GDAL_LIBRARY})

if (${BUILD_APPS})
	message(STATUS ${CGAL_LIBRARY})
	file(GLOB PC_COMPUTERS src/pc/computers/*.cpp)
	add_library (geopc SHARED src/geo.cpp src/pointcloud.cpp src/pc/normalizer.cpp src/pc/rasterizer_extsort.cpp 
		src/pc/filter.cpp ${PC_COMPUTERS})
	target_link_libraries (geopc geoutil georaster geoann ${LIBLAS_LIBRARY} CGAL CGAL_Core)

endif()

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

	# add_executable (polygonize src/apps/polygonize.cpp)
	# target_link_libraries (polygonize geoutil georaster ${Boost_LIBRARIES})

	#add_executable (pc2grid src/apps/pc2grid.cpp)
	#target_link_libraries (pc2grid geoann geoutil georaster geopc geoann ${GEOS_LIBRARY})

	# add_executable (pc2tile src/apps/pc2tile.cpp)
	# target_link_libraries (pc2tile geoutil geopc ${Boost_LIBRARIES} ${LIBLAS_LIBRARY} ${GEOS_LIBRARY})

	# add_executable (pcnorm src/apps/pcnorm.cpp)
	# target_link_libraries (pcnorm geoutil geopc ${Boost_LIBRARIES} ${LIBLAS_LIBRARY} ${GEOS_LIBRARY})

	# add_executable (pcperturb src/apps/pcperturb.cpp)
	# target_link_libraries (pcperturb geoutil geopc ${Boost_LIBRARIES} ${LIBLAS_LIBRARY})

	# add_executable (pcregister src/apps/pcregister.cpp src/pc/trajectory.cpp)
	# target_link_libraries (pcregister geoutil geopc ${Boost_LIBRARIES} ${LIBLAS_LIBRARY} ${PCL_LIBRARIES})

	# add_executable (geoinfo src/apps/geoinfo.cpp)
	# target_link_libraries (geoinfo geoutil ${GDAL_LIBRARY} ${LIBLAS_LIBRARY} ${JSONCPP_LIBRARY})

	add_executable (voidfill src/apps/voidfill.cpp)
	target_link_libraries (voidfill geoutil georaster ${GDAL_LIBRARY})

	# add_executable (datumshift src/apps/datumshift.cpp)
	# target_link_libraries (datumshift geoutil georaster ${GDAL_LIBRARY} ${PROJ_LIBRARY})

	# add_executable (extmergesort src/apps/externalmergesort.cpp)
	# target_link_libraries (extmergesort geoutil geopc ${LIBLAS_LIBRARY})

	# add_executable (sbet2csv src/apps/sbet2csv.cpp)
	# target_link_libraries (sbet2csv geoutil ${PROJ_LIBRARY})

	# add_executable (groundassign src/apps/groundassign.cpp)
	# target_link_libraries (groundassign ${LIBLAS_LIBRARY})

	add_executable (rastermatch src/apps/rastermatch.cpp)
	target_link_libraries (rastermatch geoutil geoann georaster ${GEOS_LIBRARY} ${GDAL_LIBRARY})

	#add_executable (rbf src/tests/rbf.cpp)
	#target_link_libraries (rbf)

	#install(TARGETS pc2grid pc2tile rastermatch geopc georaster geoinfo geoutil geoann RUNTIME DESTINATION bin LIBRARY DESTINATION lib) 

else()

	install(TARGETS georaster geoutil geodb RUNTIME DESTINATION bin LIBRARY DESTINATION lib) 

endif()
