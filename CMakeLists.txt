project(digital-surface-regularization-TVCG)

CMAKE_MINIMUM_REQUIRED(VERSION 3.11)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Polyscope---------------------------------------------------------------------
add_subdirectory("${PROJECT_SOURCE_DIR}/deps/polyscope")

# Eigen ------------------------------------------------------------------------
# Find Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# DGtal ------------------------------------------------------------------------
SET(BUILD_EXAMPLES OFF CACHE BOOL "DGtal examples")
add_subdirectory("${PROJECT_SOURCE_DIR}/deps/DGtal")
include("${PROJECT_BINARY_DIR}/deps/DGtal/DGtalConfig.cmake")
include_directories("${DGTAL_INCLUDE_DIRS}")

## Main target
add_executable(regularization main.cpp)
target_link_libraries(regularization polyscope DGtal ${DGTAL_LIBRARIES} Eigen3::Eigen)
target_include_directories(regularization PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/libigl/include")
