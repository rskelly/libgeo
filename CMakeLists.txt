cmake_minimum_required(VERSION 3.1)

enable_language(C)
enable_language(CXX)
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -D_GLIBCXX_PARALLEL")

if (${APPLE})
	set (GDAL_LIBRARY /usr/local/opt/gdal2/lib/libgdal.dylib)
	set (GDAL_INCLUDE_DIR /usr/local/opt/gdal2/include)
	set (CMAKE_C_COMPILER /usr/local/bin/clang-omp CACHE STRING "C compiler" FORCE)
	set (CMAKE_CXX_COMPILER /usr/local/bin/clang-omp++ CACHE STRING "C++ compiler" FORCE)
	set (CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
elseif(${WIN32})
	set (Boost_COMPILER -vc140)
	set (Boost_NO_SYSTEM_PATHS  ON)
	set (BOOST_INCLUDEDIR "E:\\boost_1_63_0")
	set (BOOST_LIBRARYDIR "E:\\boost_1_63_0\\lib64-msvc-14.0")
	set (LIB_PATH "E:\\opt\\lib")
	set (INC_PATH "E:\\opt\\include")
	set (GDAL_INCLUDE_DIR "E:\\opt\\include")
	set (GDAL_LIBRARY "E:\\opt\\lib\\gdal_i.lib")
else()
	option(CGAL_DISABLE_ROUNDING_MATH_CHECK "Disable rounding math check in CGAL. This permits Valgrind to run." ON)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
endif()

set (Boost_DEBUG OFF)
set (Boost_USE_STATIC_LIBS OFF)
set (Boost_USE_STATIC_RUNTIME OFF)

find_package (OpenMP REQUIRED)
find_package (Boost COMPONENTS system filesystem date_time REQUIRED)
find_package (GDAL 2.1 REQUIRED)
set (GEOS_LIBRARY geos)

if(OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
	set (CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
endif()


project (libgeo)

# Configure directories

link_directories(. ${Boost_LIBRARY_DIRS} ${LIB_PATH})

# Build support libraries ################################################################################

add_library (geo STATIC src/geo.cpp)
target_include_directories(geo PUBLIC include)


add_library (util STATIC src/util.cpp)
target_include_directories(util PUBLIC include ${Boost_INCLUDE_DIRS} ${GDAL_INCLUDE_DIR})
target_link_libraries (util ${Boost_LIBRARIES} ${GDAL_LIBRARY} geo rt)

add_library (db STATIC src/db.cpp)
target_include_directories(db PUBLIC include ${Boost_INCLUDE_DIRS} ${GDAL_INCLUDE_DIR})
target_link_libraries (db ${GDAL_LIBRARY} util geo)

# Build tests ############################################################################################

add_library (raster STATIC src/raster.cpp)
target_include_directories(raster PUBLIC include ${Boost_INCLUDE_DIRS} ${GDAL_INCLUDE_DIR})
target_link_libraries (raster ${GDAL_LIBRARY} ${GEOS_LIBRARY} util geo db)

add_executable (polygonize src/tests/polygonize.cpp)
target_include_directories(polygonize PUBLIC include)
target_link_libraries (polygonize PUBLIC raster)

