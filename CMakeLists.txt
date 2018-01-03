cmake_minimum_required(VERSION 3.1)

enable_language(C)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -D_GLIBCXX_PARALLEL")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")

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
	include_directories (BEFORE SYSTEM include)
	include_directories (BEFORE SYSTEM ann/include)
	include_directories (BEFORE SYSTEM /usr/include/eigen3)
	include_directories (BEFORE SYSTEM /usr/include/gdal)
	include_directories (BEFORE SYSTEM /opt/include)
endif()

set (Boost_DEBUG OFF)
set (Boost_USE_STATIC_LIBS OFF)
set (Boost_USE_STATIC_RUNTIME OFF)

find_package (OpenMP REQUIRED)
find_package (Boost COMPONENTS system filesystem date_time REQUIRED)
find_package (GDAL 2.2 REQUIRED)

find_library (GEOS_LIBRARY geos PATHS /opt/lib NO_DEFAULT_PATH)
find_library (CRYPTO_LIB crypto)

if(OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
	set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
endif()


project (libgeo)


# Configure directories###################################################################################

link_directories(. ${Boost_LIBRARY_DIRS} ${LIB_PATH})

# Build support libraries ################################################################################

add_library (geo SHARED src/geo.cpp)

add_library (geoutil SHARED src/util.cpp src/md5.cpp)
target_link_libraries (geoutil geo rt ${Boost_LIBRARIES} ${GDAL_LIBRARY} ${CRYPTO_LIB})

add_library (geodb SHARED src/db.cpp)
target_link_libraries (geodb geo geoutil ${GDAL_LIBRARY})

add_library (georaster SHARED src/raster.cpp)
target_link_libraries (georaster geo geoutil geodb ${GEOS_LIBRARY} ${GDAL_LIBRARY})

add_library (geoann SHARED ann/src/ANN.cpp ann/src/brute.cpp ann/src/kd_tree.cpp ann/src/kd_util.cpp 
	ann/src/kd_split.cpp ann/src/kd_dump.cpp ann/src/kd_search.cpp ann/src/kd_pr_search.cpp 
	ann/src/kd_fix_rad_search.cpp ann/src/bd_tree.cpp ann/src/bd_search.cpp ann/src/bd_pr_search.cpp 
	ann/src/bd_fix_rad_search.cpp ann/src/perf.cpp)

# Tests ##################################################################################################

add_executable (kdtree_test src/tests/kdtree_test.cpp)
target_link_libraries (kdtree_test geoann ${Boost_LIBRARIES})

add_executable (polygonize src/apps/polygonize.cpp)
target_link_libraries (polygonize geo geoutil georaster ${Boost_LIBRARIES})