cmake_minimum_required (VERSION 3.5)
project (lib_cpp_hybrid_a_star)

set(CMAKE_BUILD_TYPE "DEBUG")
set(CMAKE_CXX_FLAGS "-std=c++17")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g --debug-cpp --trace-expand")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_MATH_DISABLE_FLOAT128 -O1 -Wall -g --debug")


find_package(Eigen3 REQUIRED)
# find python libraries
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)
find_package(PythonLibs 3.0 REQUIRED)
include_directories(${PYTHON3_INCLUDE_DIRS} ${NumPy_INCLUDE_DIRS})

# populate matplotlib repository
include(FetchContent)
FetchContent_Declare(
    matplotlib
    GIT_REPOSITORY https://github.com/lava/matplotlib-cpp.git
    GIT_TAG        f23347fca25219d1c42cbb91608b5556814bf572
)
FetchContent_GetProperties(matplotlib)
if(NOT matplotlib_POPULATED)
    FetchContent_Populate(matplotlib)
endif()

include_directories( 
    SYSTEM 
    include
    ${Eigen3_INCLUDE_DIRS}
    ${matplotlib_SOURCE_DIR}
)

add_library(lib_cpp_hybrid_a_star src/rs_paths.cpp src/grid_a_star.cpp src/trailerlib.cpp src/trailer_hybrid_a_star.cpp)
target_link_libraries(lib_cpp_hybrid_a_star  Eigen3::Eigen ${PYTHON_LIBRARIES} Python3::NumPy)

add_executable(hybrid_a_star src/main.cpp)
target_link_libraries(hybrid_a_star  Eigen3::Eigen lib_cpp_hybrid_a_star  ${PYTHON_LIBRARIES} Python3::NumPy)
