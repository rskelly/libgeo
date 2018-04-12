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

if (${APPLE})
	set (CMAKE_C_COMPILER /usr/local/bin/clang-omp CACHE STRING "C compiler" FORCE)
	set (CMAKE_CXX_COMPILER /usr/local/bin/clang-omp++ CACHE STRING "C++ compiler" FORCE)
	set (CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
elseif(${WIN32})
	set (Boost_COMPILER -vc140)
	set (Boost_NO_SYSTEM_PATHS ON)
	set (BOOST_INCLUDEDIR "C:\\dev\\boost_1_64_0")
	set (Boost_INCLUDE_DIRS "C:\\dev\\boost_1_64_0")
	set (BOOST_LIBRARYDIR "C:\\dev\\boost_1_64_0\\lib64-msvc-14.0")
	set (CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};C:\\Qt\\5.9.1\\msvc2015_64")
	set (LIB_PATH "C:\\dev\\geos\\makefiles\\lib\\Release;C:\\dev\\gisinternals\\lib")
	set (INC_PATH "C:\\dev\\geos\\include;C:\\dev\\gisinternals\\include;C:\\dev\\eigen3")
	set (GDAL_INCLUDE_DIR "C:\\dev\\gisinternals\\include")
	set (GDAL_LIBRARY "C:\\dev\\gisinternals\\lib\\gdal_i.lib")
	set (libLAS_INCLUDE_DIR "C:\\dev\\gisinternals\\include")
	set (libLAS_LIBRARY_DIRS "C:\\dev\\gisinternals\\lib")
else()
#	include_directories (BEFORE SYSTEM include ann/include /usr/include /usr/include/eigen3 
#		/usr/include/gdal /usr/include/jsoncpp /opt/include)
endif()

set (Boost_DEBUG OFF)
set (Boost_USE_STATIC_LIBS OFF)
set (Boost_USE_STATIC_RUNTIME OFF)

set (THREADS_PREFER_PTHREAD_FLAG ON)
find_package (Threads REQUIRED)
find_package (OpenMP REQUIRED)
find_package (Boost 1.58 COMPONENTS system filesystem date_time thread REQUIRED)
find_package (GDAL 2.2 REQUIRED)
find_package (CGAL REQUIRED COMPONENTS Core)
#find_package (PCL 1.8 REQUIRED)

find_library (LIBLAS_LIBRARY las)
find_library (JSONCPP_LIBRARY jsoncpp)
find_library (CRYPTO_LIB crypto PATHS /usr/lib/x86_64-linux-gnu)
find_library (GEOS_LIBRARY geos PATHS /opt/lib)
find_library (PROJ_LIBRARY proj)

message(${GEOS_LIBRARY})

if(OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
	set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
endif()


project (libgeo)


# Configure directories###################################################################################

include_directories(BEFORE SYSTEM include ann/include /usr/include /usr/include/eigen3 
		/usr/include/gdal /usr/include/jsoncpp /opt/include ${PCL_INCLUDE_DIRS})
link_directories(BEFORE SYSTEM . /opt/lib ${Boost_LIBRARY_DIRS} ${LIB_PATH} ${PCL_LIBRARY_DIRS})

# Build support libraries ################################################################################

add_library (geoutil SHARED src/geo.cpp src/util.cpp src/md5.cpp)
target_link_libraries (geoutil rt ${Boost_LIBRARIES} ${GDAL_LIBRARY} ${CRYPTO_LIB})

add_library (geodb SHARED src/geo.cpp src/db.cpp)
target_link_libraries (geodb geoutil)

add_library (georaster SHARED src/geo.cpp src/raster.cpp src/rast/polygonize.cpp)
target_link_libraries (georaster geoutil geodb Threads::Threads ${GEOS_LIBRARY})

file(GLOB PC_COMPUTERS src/pc/computers/*.cpp)
add_library (geopc SHARED src/geo.cpp src/pointcloud.cpp src/pc/normalizer.cpp src/pc/rasterizer_extsort.cpp 
	src/pc/filter.cpp ${PC_COMPUTERS})
target_link_libraries (geopc geoutil georaster geoann ${LIBLAS_LIBRARY} ${CGAL_LIBRARY})

add_library (geoann SHARED src/geo.cpp ann/src/ANN.cpp ann/src/brute.cpp ann/src/kd_tree.cpp ann/src/kd_util.cpp 
	ann/src/kd_split.cpp ann/src/kd_dump.cpp ann/src/kd_search.cpp ann/src/kd_pr_search.cpp 
	ann/src/kd_fix_rad_search.cpp ann/src/bd_tree.cpp ann/src/bd_search.cpp ann/src/bd_pr_search.cpp 
	ann/src/bd_fix_rad_search.cpp ann/src/perf.cpp)
target_link_libraries (geoann geoutil)

# Applications ##################################################################################################

add_executable (kdtree_test src/tests/kdtree_test.cpp)
target_link_libraries (kdtree_test geoann ${Boost_LIBRARIES})

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

add_executable (sbet2csv src/apps/sbet2csv.cpp)
target_link_libraries (sbet2csv geoutil ${PROJ_LIBRARY})

add_executable (extmergesort src/apps/externalmergesort.cpp)
target_link_libraries (extmergesort geoutil geopc ${LIBLAS_LIBRARY})


install(TARGETS pc2grid pc2tile geopc georaster geoinfo geoutil geoann RUNTIME DESTINATION bin LIBRARY DESTINATION lib) 