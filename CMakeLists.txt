cmake_minimum_required(VERSION 3.10)

project(mdsim
    VERSION 1.0.0
    DESCRIPTION "Deterministic, high-performance molecular dynamics simulator"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Determinism: pin a fast-math-free, no-contraction flag set.
# -ffp-contract=off / -ffloat-store guard against FMA and excess precision,
# though the physics core uses integer fixed-point arithmetic exclusively.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-ffp-contract=off -fno-fast-math)
    if(NOT MSVC)
        add_compile_options(-ffloat-store)
    endif()
elseif(MSVC)
    add_compile_options(/fp:strict)
endif()

option(MDSIM_ENABLE_THREADS "Enable multi-threaded domain decomposition" ON)
option(MDSIM_BUILD_TESTS "Build the test suite" ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# ---------------------------------------------------------------------------
# Core library (header-only math + compiled simulation core)
# ---------------------------------------------------------------------------
add_library(mdsim_core
    src/fixed.cpp
    src/rng.cpp
    src/transcendental.cpp
    src/cell_list.cpp
    src/neighbor_list.cpp
    src/forces.cpp
    src/ewald.cpp
    src/integrator.cpp
    src/thermostat.cpp
    src/barostat.cpp
    src/shake.cpp
    src/domain.cpp
    src/thread_pool.cpp
    src/checkpoint.cpp
    src/analysis.cpp
    src/sha256.cpp
    src/particles.cpp
    src/simulation.cpp
)

target_include_directories(mdsim_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

if(MDSIM_ENABLE_THREADS)
    target_compile_definitions(mdsim_core PUBLIC MDSIM_ENABLE_THREADS=1)
endif()

# ---------------------------------------------------------------------------
# Main executable
# ---------------------------------------------------------------------------
add_executable(mdsim src/main.cpp)
target_link_libraries(mdsim PRIVATE mdsim_core)

if(MINGW)
    # Ensure the MinGW runtime DLLs are found next to the binary.
    set_target_properties(mdsim PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
    target_compile_definitions(mdsim_core PUBLIC WIN32_LEAN_AND_MEAN NOMINMAX)
endif()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(MDSIM_BUILD_TESTS)
    enable_testing()

    add_executable(test_fixed tests/test_fixed.cpp)
    target_link_libraries(test_fixed PRIVATE mdsim_core)
    add_test(NAME fixed COMMAND test_fixed)

    add_executable(test_rng tests/test_rng.cpp)
    target_link_libraries(test_rng PRIVATE mdsim_core)
    add_test(NAME rng COMMAND test_rng)

    add_executable(test_neighbor tests/test_neighbor.cpp)
    target_link_libraries(test_neighbor PRIVATE mdsim_core)
    add_test(NAME neighbor COMMAND test_neighbor)

    add_executable(test_repro tests/test_repro.cpp)
    target_link_libraries(test_repro PRIVATE mdsim_core)
    add_test(NAME repro COMMAND test_repro)
endif()

install(TARGETS mdsim RUNTIME DESTINATION bin)
