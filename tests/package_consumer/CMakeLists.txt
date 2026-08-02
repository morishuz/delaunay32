cmake_minimum_required(VERSION 3.16)
project(delaunay32_package_consumer LANGUAGES CXX)

find_package(delaunay32 CONFIG REQUIRED)

add_executable(delaunay32_package_consumer main.cpp)
target_compile_features(delaunay32_package_consumer PRIVATE cxx_std_17)
target_link_libraries(
    delaunay32_package_consumer PRIVATE delaunay32::delaunay32)

enable_testing()
add_test(NAME package_consumer COMMAND delaunay32_package_consumer)
