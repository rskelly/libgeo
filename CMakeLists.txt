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
	set (BOOST_INCLUDEDIR "C:\\boost_1_63_0")
	set (Boost_INCLUDE_DIRS "C:\\boost_1_63_0")
	set (BOOST_LIBRARYDIR "C:\\boost_1_63_0\\stage\\lib")
	set (LIB_PATH "C:\\treetops\\gisinternals\\lib")
	set (INC_PATH "C:\\treetops\\gisinternals\\include;C:\\treetops\\eigen")
	set (GDAL_INCLUDE_DIR "C:\\treetops\\gisinternals\\include")
	set (GDAL_LIBRARY "C:\\treetops\\gisinternals\\lib\\gdal_i.lib")
else()
	option(CGAL_DISABLE_ROUNDING_MATH_CHECK "Disable rounding math check in CGAL. This permits Valgrind to run." ON)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
	set (LIB_PATH /usr/lib/x86_64-linux-gnu)
	set (INC_PATH /usr/include)
endif()

# disallow in-source builds
# if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
#  message(SEND_ERROR "In-source builds are not allowed.")
# endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# add_subdirectory(stxxl)
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")

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


add_library (util STATIC src/util.cpp src/md5.cpp)
target_include_directories(util PUBLIC include ${INC_PATH} ${Boost_INCLUDE_DIRS} ${GDAL_INCLUDE_DIR})
if (UNIX AND NOT APPLE)
	target_link_libraries (util ${Boost_LIBRARIES} ${GDAL_LIBRARY} geo rt)
else ()
	target_link_libraries (util ${Boost_LIBRARIES} ${GDAL_LIBRARY} geo)
endif ()

add_library (db STATIC src/db.cpp)
target_include_directories(db PUBLIC include ${INC_PATH} ${Boost_INCLUDE_DIRS} ${GDAL_INCLUDE_DIR})
target_link_libraries (db ${GDAL_LIBRARY} util geo)

add_library (raster STATIC src/raster.cpp)
target_include_directories(raster PUBLIC include ${INC_PATH} ${Boost_INCLUDE_DIRS} ${GDAL_INCLUDE_DIR})
target_link_libraries (raster ${GDAL_LIBRARY} ${GEOS_LIBRARY} util geo db)

# Build tests ############################################################################################

add_executable (polygonize_test src/tests/polygonize.cpp)
target_include_directories(polygonize_test PUBLIC include)
target_link_libraries (polygonize_test PUBLIC raster)

add_executable (smooth_test src/tests/smooth.cpp)
target_include_directories(smooth_test PUBLIC include)
target_link_libraries (smooth_test PUBLIC raster)

add_executable (fqtree_test src/tests/fqtree.cpp)
target_include_directories(fqtree_test PUBLIC include)# ${STXXL_INCLUDE_DIRS})
target_link_libraries (fqtree_test PUBLIC util)# ${STXXL_LIBRARIES})

add_executable (mqtree_test src/tests/mqtree.cpp)
target_include_directories(mqtree_test PUBLIC include)# ${STXXL_INCLUDE_DIRS})
target_link_libraries (mqtree_test PUBLIC util)# ${STXXL_LIBRARIES})

