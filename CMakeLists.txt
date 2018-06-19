cmake_minimum_required(VERSION 3.1)

set (CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

enable_language(C)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -D_GLIBCXX_PARALLEL")

if (${APPLE})
	set (CMAKE_C_COMPILER /usr/local/bin/clang-omp CACHE STRING "C compiler" FORCE)
	set (CMAKE_CXX_COMPILER /usr/local/bin/clang-omp++ CACHE STRING "C++ compiler" FORCE)
	set (CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
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
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
	set (GEOS_INCLUDEDIR /opt/include)
	set (GEOS_LIBRARYDIR /opt/lib)
	set (EIGEN_INCLUDEDIR /usr/include/eigen3)
endif()

set (Boost_DEBUG OFF)
set (Boost_USE_SHARED_LIBS OFF)
set (Boost_USE_SHARED_RUNTIME OFF)

find_package (OpenMP REQUIRED)
find_package (Boost COMPONENTS system filesystem date_time REQUIRED)
find_package (GDAL 2.1 REQUIRED)

find_library (GEOS_LIBRARY geos PATHS ${GEOS_LIBRARYDIR})
message (STATUS "GEOS: ${GEOS_LIBRARY} ${GEOS_LIBRARYDIR}")

if(OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
	set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
endif()


project (libgeo)


# Configure directories

include_directories (include ${GEOS_INCLUDEDIR} ${EIGEN_INCLUDEDIR} ${Boost_INCLUDE_DIRS} ${GDAL_INCLUDE_DIR})
link_directories(. ${GEOS_LIBRARYDIR} ${Boost_LIBRARY_DIRS} ${LIB_PATH})

# Build support libraries ################################################################################

add_library (util SHARED src/util.cpp src/md5.cpp src/geo.cpp)
if (UNIX AND NOT APPLE)
	target_link_libraries (util ${Boost_LIBRARIES} ${GDAL_LIBRARY} rt)
else ()
	target_link_libraries (util ${Boost_LIBRARIES} ${GDAL_LIBRARY})
endif ()

add_library (db SHARED src/db.cpp)
target_link_libraries (db ${GDAL_LIBRARY} util)

add_library (raster SHARED src/raster.cpp)
target_link_libraries (raster ${GDAL_LIBRARY} ${GEOS_LIBRARY} util db)

# Build tests ############################################################################################

#add_executable (polygonize src/apps/polygonize.cpp)
#target_include_directories(polygonize PUBLIC include)
#target_link_libraries (polygonize PUBLIC raster)

#add_executable (smooth_test src/tests/smooth.cpp)
#target_include_directories(smooth_test PUBLIC include)
#target_link_libraries (smooth_test PUBLIC raster)

#add_executable (fqtree_test src/tests/fqtree.cpp)
#target_include_directories(fqtree_test PUBLIC include)# ${STXXL_INCLUDE_DIRS})
#target_link_libraries (fqtree_test PUBLIC util)# ${STXXL_LIBRARIES})

#add_executable (mqtree_test src/tests/mqtree.cpp)
#target_include_directories(mqtree_test PUBLIC include)# ${STXXL_INCLUDE_DIRS})
#target_link_libraries (mqtree_test PUBLIC util)# ${STXXL_LIBRARIES})

