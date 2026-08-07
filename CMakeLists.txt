cmake_minimum_required(VERSION 3.16)
project(delaunay32 VERSION 0.6.2 LANGUAGES CXX)

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
    "Build the standalone performance benchmark"
    ${DELAUNAY32_TOP_LEVEL})
option(
    DELAUNAY32_BUILD_EXTRAS
    "Build the optional sampling, geometry JSON, and SVG companion library"
    ON)
option(
    DELAUNAY32_BUILD_EXAMPLES
    "Build the example programs"
    ${DELAUNAY32_TOP_LEVEL})
option(
    DELAUNAY32_BUILD_TESTS
    "Build the validation suite"
    ${DELAUNAY32_TOP_LEVEL})

if(DELAUNAY32_TOP_LEVEL AND
   NOT CMAKE_CONFIGURATION_TYPES AND
   NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

function(delaunay32_enable_warnings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4)
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        target_compile_options(
            ${target}
            PRIVATE
                -Wall
                -Wextra
                -Wconversion
                -Wshadow
                -Wsign-conversion
        )
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra)
    endif()
endfunction()

set(DELAUNAY32_LIBRARY_SOURCES
    src/constraints.cpp
    src/delaunay.cpp
    src/export.cpp
    src/internal.hpp
    src/polygon.cpp
    src/quantization.cpp
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

if(DELAUNAY32_BUILD_EXTRAS)
    add_library(
        delaunay32_extras
        extras/geometry.cpp
        extras/internal.hpp
        extras/json.cpp
        extras/sampling.cpp
        extras/svg.cpp
    )
    add_library(delaunay32::extras ALIAS delaunay32_extras)
    set_target_properties(delaunay32_extras PROPERTIES EXPORT_NAME extras)
    target_compile_features(delaunay32_extras PUBLIC cxx_std_17)
    target_include_directories(delaunay32_extras
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
    target_link_libraries(
        delaunay32_extras PUBLIC delaunay32::delaunay32)
    delaunay32_enable_warnings(delaunay32_extras)
endif()

if(DELAUNAY32_BUILD_EXAMPLES)
    if(NOT DELAUNAY32_BUILD_EXTRAS)
        message(FATAL_ERROR
            "DELAUNAY32_BUILD_EXAMPLES requires "
            "DELAUNAY32_BUILD_EXTRAS=ON")
    endif()
    add_executable(
        delaunay32_logo_polygon_example
        examples/delaunay32_logo_polygon_example.cpp)
    target_link_libraries(
        delaunay32_logo_polygon_example PRIVATE delaunay32::extras)
    target_compile_definitions(
        delaunay32_logo_polygon_example
        PRIVATE
            DELAUNAY32_LOGO_INPUT_PATH="${CMAKE_CURRENT_SOURCE_DIR}/examples/data/delaunay32_logo.json")
    delaunay32_enable_warnings(delaunay32_logo_polygon_example)

    add_executable(
        delaunay32_hello_mesh_example
        examples/delaunay32_hello_mesh_example.cpp)
    target_link_libraries(
        delaunay32_hello_mesh_example PRIVATE delaunay32::extras)
    delaunay32_enable_warnings(delaunay32_hello_mesh_example)

    add_executable(
        delaunay32_quantization_example
        examples/delaunay32_quantization_example.cpp)
    target_link_libraries(
        delaunay32_quantization_example PRIVATE delaunay32::extras)
    delaunay32_enable_warnings(delaunay32_quantization_example)

    add_executable(
        delaunay32_constraints_example
        examples/delaunay32_constraints_example.cpp)
    target_link_libraries(
        delaunay32_constraints_example PRIVATE delaunay32::extras)
    delaunay32_enable_warnings(delaunay32_constraints_example)

    add_executable(
        delaunay32_polygon_example
        examples/delaunay32_polygon_example.cpp)
    target_link_libraries(
        delaunay32_polygon_example PRIVATE delaunay32::extras)
    delaunay32_enable_warnings(delaunay32_polygon_example)
endif()

if(DELAUNAY32_BUILD_BENCHMARKS)
    add_executable(delaunay_benchmark benchmarks/benchmark.cpp)
    target_link_libraries(
        delaunay_benchmark
        PRIVATE delaunay32::delaunay32)
    delaunay32_enable_warnings(delaunay_benchmark)
endif()

if(DELAUNAY32_BUILD_TESTS)
    enable_testing()
    add_executable(delaunay_tests tests/validation.cpp)
    target_include_directories(
        delaunay_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks)
    target_link_libraries(
        delaunay_tests
        PRIVATE delaunay32::delaunay32)
    delaunay32_enable_warnings(delaunay_tests)
    add_test(NAME validation COMMAND delaunay_tests)

    if(DELAUNAY32_BUILD_EXTRAS)
        add_executable(
            delaunay_extras_tests
            tests/extras_validation.cpp)
        target_link_libraries(
            delaunay_extras_tests PRIVATE delaunay32::extras)
        delaunay32_enable_warnings(delaunay_extras_tests)
        add_test(
            NAME extras_validation
            COMMAND delaunay_extras_tests
                    "${CMAKE_CURRENT_BINARY_DIR}/extras-roundtrip.json"
                    "${CMAKE_CURRENT_BINARY_DIR}/extras-mesh.svg"
                    "${CMAKE_CURRENT_BINARY_DIR}/extras-polygon.svg")
    endif()

    if(DELAUNAY32_BUILD_EXAMPLES)
        add_test(
            NAME example_hello_mesh
            COMMAND delaunay32_hello_mesh_example)
        add_test(
            NAME example_quantization
            COMMAND delaunay32_quantization_example)
        add_test(
            NAME example_constraints
            COMMAND delaunay32_constraints_example)
        add_test(
            NAME example_polygon
            COMMAND delaunay32_polygon_example)
        add_test(
            NAME example_logo_polygon
            COMMAND
                delaunay32_logo_polygon_example
                "${CMAKE_CURRENT_SOURCE_DIR}/examples/data/delaunay32_logo.json"
                "${CMAKE_CURRENT_BINARY_DIR}/example-logo-polygon.svg")
        set_tests_properties(
            example_hello_mesh
            example_quantization
            example_constraints
            example_polygon
            example_logo_polygon
            PROPERTIES LABELS examples)
    endif()

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

endif()

set(DELAUNAY32_INSTALL_TARGETS delaunay32)
if(DELAUNAY32_BUILD_EXTRAS)
    list(APPEND DELAUNAY32_INSTALL_TARGETS delaunay32_extras)
endif()
install(
    TARGETS ${DELAUNAY32_INSTALL_TARGETS}
    EXPORT delaunay32Targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(
    FILES
        include/delaunay32/delaunay.hpp
        include/delaunay32/quantization.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/delaunay32
)
if(DELAUNAY32_BUILD_EXTRAS)
    install(
        DIRECTORY include/delaunay32/extras
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/delaunay32
    )
endif()
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
