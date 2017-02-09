cmake_minimum_required(VERSION 3.1)

set (CMAKE_POSITION_INDEPENDENT_CODE ON)

set (Boost_DEBUG OFF)
set (Boost_USE_STATIC_LIBS OFF)
set (Boost_USE_STATIC_RUNTIME OFF)

set (GEOS_LIBRARY geos)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -D_GLIBCXX_PARALLEL")

if (${APPLE})

	# For GDAL2
	set (GDAL_LIBRARY /usr/local/opt/gdal2/lib/libgdal.dylib)
	set (GDAL_INCLUDE_DIR /usr/local/opt/gdal2/include)
	set (CMAKE_C_COMPILER /usr/local/bin/clang-omp CACHE STRING "C compiler" FORCE)
	set (CMAKE_CXX_COMPILER /usr/local/bin/clang-omp++ CACHE STRING "C++ compiler" FORCE)
	set (CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")

elseif(${WIN32})
	
	# Preparing to build on Windows
	# 1) download and install vs 2015
	# 2) use web installer to install qt5
	# 3) download prebuild boost for vs msvc-14 -- build using bootstrap, b2
	# 3) download gisinternals lib and headers zips for msvc-12 -- unzip to C:\opt
	# 4) download and install cmake
	# 5) download and install git
	# --- 6) checkout liblas from git
	# 7) download eigen -- extract to c:\opt\include -- rename folder to eigen3
	
	# Configure boost.
	set (Boost_COMPILER -vc140)
	set (Boost_NO_SYSTEM_PATHS  ON)
	set (BOOST_INCLUDEDIR "E:\\boost_1_63_0")
	set (BOOST_LIBRARYDIR "E:\\boost_1_63_0\\lib64-msvc-14.0")
	set (CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};C:\\Qt\\5.7\\msvc2015_64\\lib\\cmake\\Qt5Widgets")
	set (LIB_PATH "E:\\opt\\lib")
	set (INC_PATH "E:\\opt\\include")
	set (GDAL_INCLUDE_DIR "E:\\opt\\include")
	set (GDAL_LIBRARY "E:\\opt\\lib\\gdal_i.lib")
	set (libLAS_INCLUDE_DIR "E:\\opt\\include")
	set (libLAS_LIBRARY_DIRS "E:\\opt\\lib")
	#set (libLAS_BINARY_DIRS "E:\\opt\\bin")
	
else()
	option(CGAL_DISABLE_ROUNDING_MATH_CHECK "Disable rounding math check in CGAL. This permits Valgrind to run." ON)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
endif()

enable_language(C)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

set (CMAKE_AUTOUIC ON)
set (CMAKE_AUTOMOC ON)
set (CMAKE_INCLUDE_CURRENT_DIR ON)
set (CMAKE_AUTOMOC_OPTIONS "-Iinclude")

find_package (OpenMP REQUIRED)
find_package (Boost COMPONENTS system filesystem date_time REQUIRED)
find_package (GDAL 2.1 REQUIRED)
find_package (libLAS 1.8) # Will fail on windows with gisinternals linkage; ignore.

if(OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
	set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

project (geolibs)

# Configure directories
link_directories(. ${Boost_LIBRARY_DIRS} ${LIB_PATH})
include_directories(include ${Boost_INCLUDE_DIRS} ${GDAL_INCLUDE_DIR} ${Qt5_DIR}/include ${Qt5Core_INCLUDE_DIRS} ${Qt5Widgets_INCLUDE_DIRS} ${Qt5Gui_INCLUDE_DIRS} ${INC_PATH})


# Build support libraries ############################################################################################

add_library (geotools STATIC src/geotools.cpp)

add_library (util STATIC src/util.cpp)
target_link_libraries (util ${Boost_LIBRARIES} geotools)

add_library (db STATIC src/db.cpp)
target_link_libraries (db ${GDAL_LIBRARY} geotools)

add_library (raster STATIC src/raster.cpp)
target_link_libraries (raster ${GDAL_LIBRARY} ${GEOS_LIBRARY} util geotools db)

# Install ###########################################################################################################

#