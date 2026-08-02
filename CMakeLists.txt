cmake_minimum_required(VERSION 3.16)
project(delaunay32 VERSION 0.3.0 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
find_package(Threads REQUIRED)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    set(DELAUNAY32_TOP_LEVEL ON)
else()
    set(DELAUNAY32_TOP_LEVEL OFF)
endif()

option(
    DELAUNAY32_BUILD_BENCHMARKS
    "Build the validated Delaunator comparison"
    ${DELAUNAY32_TOP_LEVEL})
option(
    DELAUNAY32_BUILD_EXAMPLES
    "Build the dependency-free SVG examples"
    ${DELAUNAY32_TOP_LEVEL})
option(
    DELAUNAY32_BUILD_TESTS
    "Build the validation suite"
    ${DELAUNAY32_TOP_LEVEL})

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

function(delaunay32_enable_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4)
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra)
    endif()
endfunction()

set(DELAUNAY32_LIBRARY_SOURCES
    src/delaunay.cpp
    src/export.cpp
    src/internal.hpp
    src/result.cpp
    src/topology.cpp
    src/worker_team.cpp
)
add_library(delaunay32 ${DELAUNAY32_LIBRARY_SOURCES})
add_library(delaunay32::delaunay32 ALIAS delaunay32)
target_compile_features(delaunay32 PUBLIC cxx_std_17)
target_include_directories(delaunay32
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(delaunay32 PUBLIC Threads::Threads)
delaunay32_enable_warnings(delaunay32)

if(DELAUNAY32_BUILD_EXAMPLES)
    add_executable(
        delaunay_svg_example
        examples/delaunay_svg_example.cpp
        examples/svg_io.cpp
        examples/svg_io.hpp
    )
    target_link_libraries(
        delaunay_svg_example PRIVATE delaunay32::delaunay32)
    delaunay32_enable_warnings(delaunay_svg_example)

    add_executable(
        delaunay_float_example
        examples/delaunay_float_example.cpp
        examples/svg_io.cpp
        examples/svg_io.hpp
    )
    target_link_libraries(
        delaunay_float_example PRIVATE delaunay32::delaunay32)
    delaunay32_enable_warnings(delaunay_float_example)
endif()

set(DELAUNATOR_HEADER
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party/delaunator-cpp/include/delaunator.hpp")
if(DELAUNAY32_BUILD_BENCHMARKS OR DELAUNAY32_BUILD_TESTS)
    if(NOT EXISTS "${DELAUNATOR_HEADER}")
        message(FATAL_ERROR
            "Delaunator submodule is missing. Run: "
            "git submodule update --init --recursive")
    endif()

    add_library(delaunator_baseline STATIC
        benchmarks/delaunator_adapter.cpp
        benchmarks/delaunator_adapter.hpp
    )
    target_compile_features(delaunator_baseline PUBLIC cxx_std_17)
    target_include_directories(delaunator_baseline
        PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks
            ${CMAKE_CURRENT_SOURCE_DIR}/include
    )
    target_include_directories(
        delaunator_baseline
        SYSTEM PUBLIC
            ${CMAKE_CURRENT_SOURCE_DIR}/third_party/delaunator-cpp/include)
    delaunay32_enable_warnings(delaunator_baseline)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        target_compile_options(
            delaunator_baseline PRIVATE -Wno-unused-but-set-variable)
    endif()
endif()

if(DELAUNAY32_BUILD_BENCHMARKS)
    add_executable(delaunay_benchmark benchmarks/benchmark.cpp)
    target_link_libraries(
        delaunay_benchmark
        PRIVATE delaunay32::delaunay32 delaunator_baseline)
    delaunay32_enable_warnings(delaunay_benchmark)
endif()

if(DELAUNAY32_BUILD_TESTS)
    enable_testing()
    add_executable(delaunay_tests tests/validation.cpp)
    target_link_libraries(
        delaunay_tests
        PRIVATE delaunay32::delaunay32 delaunator_baseline)
    delaunay32_enable_warnings(delaunay_tests)
    add_test(NAME validation COMMAND delaunay_tests)

    add_executable(
        delaunay_stress_tests
        tests/stress_validation.cpp)
    target_include_directories(
        delaunay_stress_tests
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks)
    target_link_libraries(
        delaunay_stress_tests PRIVATE delaunay32::delaunay32)
    delaunay32_enable_warnings(delaunay_stress_tests)
    add_test(NAME stress_validation COMMAND delaunay_stress_tests)
    set_tests_properties(stress_validation PROPERTIES TIMEOUT 300)

    # Compile a private variant that deterministically exhausts the fixed
    # parallel edge arena. This verifies the serial recovery path without
    # exposing test controls in the public API or relying on pathological data.
    add_library(
        delaunay32_limited_parallel_arena STATIC
        ${DELAUNAY32_LIBRARY_SOURCES})
    target_compile_features(
        delaunay32_limited_parallel_arena PUBLIC cxx_std_17)
    target_include_directories(
        delaunay32_limited_parallel_arena
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(
        delaunay32_limited_parallel_arena PUBLIC Threads::Threads)
    target_compile_definitions(
        delaunay32_limited_parallel_arena
        PRIVATE DELAUNAY32_TEST_PARALLEL_EDGE_ARENA_DART_LIMIT=262144)
    delaunay32_enable_warnings(delaunay32_limited_parallel_arena)

    add_executable(
        delaunay_parallel_fallback_tests
        tests/parallel_fallback.cpp)
    target_include_directories(
        delaunay_parallel_fallback_tests
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks)
    target_link_libraries(
        delaunay_parallel_fallback_tests
        PRIVATE delaunay32_limited_parallel_arena)
    delaunay32_enable_warnings(delaunay_parallel_fallback_tests)
    add_test(
        NAME parallel_arena_fallback
        COMMAND delaunay_parallel_fallback_tests)
    set_tests_properties(parallel_arena_fallback PROPERTIES TIMEOUT 300)

    if(DELAUNAY32_BUILD_EXAMPLES)
        add_test(
            NAME svg_example
            COMMAND delaunay_svg_example 100
                    "${CMAKE_CURRENT_BINARY_DIR}/mesh-smoke.svg" 42)
        add_test(
            NAME svg_example_csv
            COMMAND delaunay_svg_example
                    --input
                    "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/svg_points.csv"
                    --output
                    "${CMAKE_CURRENT_BINARY_DIR}/mesh-csv-smoke.svg")
        add_test(
            NAME svg_example_csv_collinear
            COMMAND delaunay_svg_example
                    --input
                    "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/svg_collinear.csv"
                    --output
                    "${CMAKE_CURRENT_BINARY_DIR}/mesh-collinear-smoke.svg")
        add_test(
            NAME float_svg_example
            COMMAND delaunay_float_example
                    "${CMAKE_CURRENT_BINARY_DIR}/mesh-float-smoke.svg")
    endif()
endif()

install(
    TARGETS delaunay32
    EXPORT delaunay32Targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(
    DIRECTORY include/delaunay32
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(
    EXPORT delaunay32Targets
    FILE delaunay32Targets.cmake
    NAMESPACE delaunay32::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/delaunay32
)

configure_package_config_file(
    cmake/delaunay32Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/delaunay32Config.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/delaunay32
)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/delaunay32ConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMinorVersion
)
install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/delaunay32Config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/delaunay32ConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/delaunay32
)
install(
    FILES LICENSE
    DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/delaunay32
)
