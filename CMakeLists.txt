cmake_minimum_required(VERSION 3.1)

enable_language(C)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas -Wall -D_GLIBCXX_PARALLEL")

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set (CMAKE_INSTALL_RPATH "/opt/lib;/usr/local/lib")
set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set (CMAKE_CXX_COMPILER_LAUNCHER ccache)

set (Boost_DEBUG OFF)
set (Boost_USE_STATIC_LIBS OFF)
set (Boost_USE_STATIC_RUNTIME OFF)
set (Boost_NO_SYSTEM_PATHS ON)
set (Boost_USE_MULTITHREADED ON)

if (${APPLE})
	set (CMAKE_C_COMPILER /usr/local/bin/clang-omp CACHE STRING "C compiler" FORCE)
	set (CMAKE_CXX_COMPILER /usr/local/bin/clang-omp++ CACHE STRING "C++ compiler" FORCE)
	set (CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
elseif(${WIN32})
	# Boost is installed from source using bjam, etc.
	set (Boost_COMPILER -vc141)
	set (BOOST_INCLUDEDIR "C:/Users/IEUser/dev/boost_1_66_0")
	set (BOOST_LIBRARYDIR "C:/Users/IEUser/dev/boost_1_66_0/stage/lib")

	# OSGeo is installed to the usual location.
	set (OSGEO_LIBRARYDIR "C:/OSGeo4W64/bin")
	set (OSGEO_INCLUDEDIR "C:/OSGeo4W64/include")
	# GDAL is in the OSGeo folder, but the latest version is
	# off in its own folder, "apps."
	set (GDAL_INCLUDE_DIR "C:/OSGeo4W64/apps/gdal-dev/include")
	set (GDAL_LIBRARYDIR "C:/OSGeo4W64/apps/gdal-dev/lib")
	set (GDAL_LIBRARY gdal_i)
	# Eigen is in the OSGeo install.
	set (EIGEN_INCLUDEDIR "${OSGEO_INCLUDEDIR}/eigen3")

	# Qt is installed minimally in the usual location.
	set (QT_LIBRARYDIR "C:/Qt/5.10.1/msvc2017_64")
	
	# GEOS is built and installed from source, must be version
	# 3.6.x, but this is likely to change.
	set (GEOS_INCLUDEDIR "C:/Users/IEUser/dev/build/include")
	set (GEOS_LIBRARYDIR "C:/Users/IEUser/dev/build/lib")
	
	# Crypto is really openssl, built and installed from 
	# github.
	set (CRYPTO_LIBRARYDIR "C:/Users/IEUser/dev/git/openssl")
	set (CRYPTO_INCLUDEDIR "C:/Users/IEUser/dev/git/openssl/include")
	
	# Liblas?
	#set (libLAS_INCLUDE_DIR "C:\\dev\\gisinternals\\include")
	#set (libLAS_LIBRARY_DIRS "C:\\dev\\gisinternals\\lib")
	
	# Qt should be added to the prefix dir.
	set (CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${QT_LIBRARYDIR}")
else()
	set (EIGEN_INCLUDEDIR "/usr/include/eigen3")
	set (BOOST_INCLUDEDIR "/usr/include")
	set (BOOST_LIBRARYDIR "/usr/lib/x86_64-linux-gnu")
	set (JSON_INCLUDEDIR "/usr/include/jsoncpp")
endif()

set (THREADS_PREFER_PTHREAD_FLAG ON)

find_package (Threads REQUIRED)
find_package (OpenMP REQUIRED)
find_package (Boost 1.58 COMPONENTS system filesystem date_time thread REQUIRED)
find_package (GDAL 2.2 REQUIRED)
find_package (CGAL REQUIRED COMPONENTS Core)
find_package (PCL 1.8 REQUIRED)

find_library (LIBLAS_LIBRARY las)
find_library (JSONCPP_LIBRARY jsoncpp)
find_library (CRYPTO_LIBRARY crypto libcrypto PATHS ${CRYPTO_LIBRARYDIR})
find_library (GEOS_LIBRARY geos PATHS ${GEOS_LIBRARYDIR})
find_library (PROJ_LIBRARY proj)

if(OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
	set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
endif()

project (libgeo)

# Configure directories###################################################################################

include_directories(BEFORE SYSTEM ./include ./ann/include /opt/include
	${GEOS_INCLUDEDIR} ${OSGEO_INCLUDEDIR} ${GDAL_INCLUDE_DIR} ${PCL_INCLUDE_DIRS} 
	${Boost_INCLUDE_DIRS} ${CRYPTO_INCLUDEDIR} ${EIGEN_INCLUDEDIR} ${JSON_INCLUDEDIR}) 
link_directories(. /opt/lib ${Boost_LIBRARY_DIRS} ${GDAL_LIBRARYDIR} ${PCL_LIBRARY_DIRS})

# Build support libraries ################################################################################

add_library (geoutil SHARED src/geo.cpp src/util.cpp src/md5.cpp)
target_link_libraries (geoutil ${Boost_LIBRARIES} ${GDAL_LIBRARY} ${CRYPTO_LIBRARY})
# is rt needed here doesn't work on windows.

add_library (geodb SHARED src/geo.cpp src/db.cpp)
target_link_libraries (geodb geoutil)

add_library (georaster SHARED src/geo.cpp src/raster.cpp src/rast/polygonize.cpp)
target_link_libraries (georaster geoutil geodb Threads::Threads ${GEOS_LIBRARY})

if (${BUILD_APPS})

	file(GLOB PC_COMPUTERS src/pc/computers/*.cpp)
	add_library (geopc SHARED src/geo.cpp src/pointcloud.cpp src/pc/normalizer.cpp src/pc/rasterizer_extsort.cpp 
		src/pc/filter.cpp ${PC_COMPUTERS})
	target_link_libraries (geopc geoutil georaster geoann ${LIBLAS_LIBRARY} ${CGAL_LIBRARY})

	add_library (geoann SHARED src/geo.cpp ann/src/ANN.cpp ann/src/brute.cpp 
		ann/src/kd_tree.cpp ann/src/kd_util.cpp 
		ann/src/kd_split.cpp ann/src/kd_dump.cpp ann/src/kd_search.cpp ann/src/kd_pr_search.cpp 
		ann/src/kd_fix_rad_search.cpp ann/src/bd_tree.cpp ann/src/bd_search.cpp ann/src/bd_pr_search.cpp 
		ann/src/bd_fix_rad_search.cpp ann/src/perf.cpp)
	target_link_libraries (geoann geoutil)

endif()

# Applications ##################################################################################################

if(${BUILD_APPS})

	add_executable (polygonize src/apps/polygonize.cpp)
	target_link_libraries (polygonize geoutil georaster ${Boost_LIBRARIES})

	add_executable (pc2grid src/apps/pc2grid.cpp)
	target_link_libraries (pc2grid geoann geoutil georaster geopc geoann ${GEOS_LIBRARY})

	add_executable (pc2tile src/apps/pc2tile.cpp)
	target_link_libraries (pc2tile geoutil geopc ${Boost_LIBRARIES} ${LIBLAS_LIBRARY} ${GEOS_LIBRARY})

	add_executable (pcnorm src/apps/pcnorm.cpp)
	target_link_libraries (pcnorm geoutil geopc ${Boost_LIBRARIES} ${LIBLAS_LIBRARY} ${GEOS_LIBRARY})

	#add_executable (pcperturb src/apps/pcperturb.cpp)
	#target_link_libraries (pcperturb geoutil geopc ${Boost_LIBRARIES} ${LIBLAS_LIBRARY})

	#add_executable (pcregister src/apps/pcregister.cpp src/pc/trajectory.cpp)
	#target_link_libraries (pcregister geoutil geopc ${Boost_LIBRARIES} ${LIBLAS_LIBRARY} ${PCL_LIBRARIES})

	add_executable (geoinfo src/apps/geoinfo.cpp)
	target_link_libraries (geoinfo geoutil ${GDAL_LIBRARY} ${LIBLAS_LIBRARY} ${JSONCPP_LIBRARY})

	add_executable (voidfill src/apps/voidfill.cpp)
	target_link_libraries (voidfill geoutil georaster ${GDAL_LIBRARY})

	add_executable (datumshift src/apps/datumshift.cpp)
	target_link_libraries (datumshift geoutil georaster ${GDAL_LIBRARY} ${PROJ_LIBRARY})

	add_executable (extmergesort src/apps/externalmergesort.cpp)
	target_link_libraries (extmergesort geoutil geopc ${LIBLAS_LIBRARY})

	add_executable (sbet2csv src/apps/sbet2csv.cpp)
	target_link_libraries (sbet2csv geoutil ${PROJ_LIBRARY})

	add_executable (groundassign src/apps/groundassign.cpp)
	target_link_libraries (groundassign ${LIBLAS_LIBRARY})

	add_executable (rastermatch src/apps/rastermatch.cpp)
	target_link_libraries (rastermatch geoutil geoann georaster ${GEOS_LIBRARY} ${GDAL_LIBRARY} ${PCL_LIBRARIES})

	#add_executable (rbf src/tests/rbf.cpp)
	#target_link_libraries (rbf)

	install(TARGETS pc2grid pc2tile rastermatch geopc georaster geoinfo geoutil geoann RUNTIME DESTINATION bin LIBRARY DESTINATION lib) 

else()

	install(TARGETS georaster geoutil geodb RUNTIME DESTINATION bin LIBRARY DESTINATION lib) 

endif()
