cmake_minimum_required(VERSION 3.16)
project(DMTransport VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(transport_grid Code/src/PhaseSpaceGrid.cpp)
target_include_directories(transport_grid PUBLIC Code/include)
target_compile_features(transport_grid PUBLIC cxx_std_14)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(transport_grid PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_library(transport_projector Code/src/PhaseSpaceProjector.cpp)
target_include_directories(transport_projector PUBLIC Code/include)
target_compile_features(transport_projector PUBLIC cxx_std_14)
target_link_libraries(transport_projector PUBLIC transport_grid)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(transport_projector PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_library(transport_grid_schema Code/src/GridSchema.cpp)
target_include_directories(transport_grid_schema PUBLIC Code/include)
target_compile_features(transport_grid_schema PUBLIC cxx_std_14)
target_link_libraries(transport_grid_schema PUBLIC transport_grid)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(transport_grid_schema PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(t04_grid_schema_probe Code/tools/t04_grid_schema_probe.cpp)
target_link_libraries(t04_grid_schema_probe PRIVATE transport_grid_schema)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(t04_grid_schema_probe PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_library(transport_reference_physics
  Code/src/physics/ScatteringPhysics.cpp
  Code/src/physics/SolarBackground.cpp)
target_include_directories(transport_reference_physics PUBLIC Code/include)
target_compile_features(transport_reference_physics PUBLIC cxx_std_14)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(transport_reference_physics PRIVATE
    -Wall -Wextra -Wpedantic)
endif()

add_library(transport_cell_geometry Code/src/PhaseSpaceCellGeometry.cpp)
target_include_directories(transport_cell_geometry PUBLIC Code/include)
target_compile_features(transport_cell_geometry PUBLIC cxx_std_14)
target_link_libraries(transport_cell_geometry PUBLIC
  transport_grid transport_reference_physics)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(transport_cell_geometry PRIVATE
    -Wall -Wextra -Wpedantic)
endif()

add_library(transport_collision_kernel Code/src/CollisionKernel.cpp)
target_include_directories(transport_collision_kernel PUBLIC Code/include)
target_compile_features(transport_collision_kernel PUBLIC cxx_std_14)
target_link_libraries(transport_collision_kernel PUBLIC
  transport_cell_geometry transport_reference_physics)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(transport_collision_kernel PRIVATE
    -Wall -Wextra -Wpedantic)
endif()

add_executable(sd_rate_diagnostic Code/tools/sd_rate_diagnostic.cpp)
target_link_libraries(sd_rate_diagnostic PRIVATE transport_reference_physics)

add_executable(t03_validation_sampler
  Code/tools/t03_validation_sampler.cpp)
target_link_libraries(t03_validation_sampler
  PRIVATE transport_reference_physics)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(t03_validation_sampler PRIVATE
    -Wall -Wextra -Wpedantic)
endif()

include(CTest)
if(BUILD_TESTING)
  add_executable(test_phase_space_grid Code/tests/test_phase_space_grid.cpp)
  target_link_libraries(test_phase_space_grid PRIVATE transport_grid)
  add_test(NAME phase_space_grid COMMAND test_phase_space_grid)
  set_tests_properties(phase_space_grid PROPERTIES TIMEOUT 10)

  add_executable(test_phase_space_projector
    Code/tests/test_phase_space_projector.cpp)
  target_link_libraries(test_phase_space_projector PRIVATE transport_projector)
  add_test(NAME phase_space_projector COMMAND test_phase_space_projector)
  set_tests_properties(phase_space_projector PROPERTIES TIMEOUT 10)

  add_executable(test_phase_space_cell_geometry
    Code/tests/test_phase_space_cell_geometry.cpp)
  target_link_libraries(test_phase_space_cell_geometry
    PRIVATE transport_cell_geometry)
  add_test(NAME phase_space_cell_geometry
    COMMAND test_phase_space_cell_geometry
            ${CMAKE_CURRENT_SOURCE_DIR}/Code/data/solar/model_agss09.dat
            ${CMAKE_CURRENT_SOURCE_DIR}/Code/data/solar/Nuclear_Data.txt)
  set_tests_properties(phase_space_cell_geometry PROPERTIES TIMEOUT 10)

  add_executable(test_collision_kernel Code/tests/test_collision_kernel.cpp)
  target_link_libraries(test_collision_kernel PRIVATE
    transport_collision_kernel)
  add_test(NAME collision_kernel
    COMMAND test_collision_kernel
            ${CMAKE_CURRENT_SOURCE_DIR}/Code/data/solar/model_agss09.dat
            ${CMAKE_CURRENT_SOURCE_DIR}/Code/data/solar/Nuclear_Data.txt)
  set_tests_properties(collision_kernel PROPERTIES TIMEOUT 10)

  add_executable(test_scattering_physics
    Code/tests/test_scattering_physics.cpp)
  target_link_libraries(test_scattering_physics
    PRIVATE transport_reference_physics)
  add_test(NAME scattering_physics COMMAND test_scattering_physics)
  set_tests_properties(scattering_physics PROPERTIES TIMEOUT 10)

  add_executable(test_solar_background
    Code/tests/test_solar_background.cpp)
  target_link_libraries(test_solar_background
    PRIVATE transport_reference_physics)
  add_test(NAME solar_background
    COMMAND test_solar_background
            ${CMAKE_CURRENT_SOURCE_DIR}/Code/data/solar/model_agss09.dat
            ${CMAKE_CURRENT_SOURCE_DIR}/Code/data/solar/Nuclear_Data.txt
            ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(solar_background PROPERTIES TIMEOUT 10)

  add_executable(test_scattering_rate
    Code/tests/test_scattering_rate.cpp)
  target_link_libraries(test_scattering_rate
    PRIVATE transport_reference_physics)
  add_test(NAME scattering_rate
    COMMAND test_scattering_rate
            ${CMAKE_CURRENT_SOURCE_DIR}/Code/data/solar/model_agss09.dat
            ${CMAKE_CURRENT_SOURCE_DIR}/Code/data/solar/Nuclear_Data.txt)
  set_tests_properties(scattering_rate PROPERTIES TIMEOUT 10)

  add_executable(test_target_sampling
    Code/tests/test_target_sampling.cpp)
  target_link_libraries(test_target_sampling
    PRIVATE transport_reference_physics)
  add_test(NAME target_sampling
    COMMAND test_target_sampling
            ${CMAKE_CURRENT_SOURCE_DIR}/Code/data/solar/model_agss09.dat
            ${CMAKE_CURRENT_SOURCE_DIR}/Code/data/solar/Nuclear_Data.txt)
  set_tests_properties(target_sampling PROPERTIES TIMEOUT 30)

  add_executable(test_single_collision
    Code/tests/test_single_collision.cpp)
  target_link_libraries(test_single_collision
    PRIVATE transport_reference_physics)
  add_test(NAME single_collision
    COMMAND test_single_collision
            ${CMAKE_CURRENT_SOURCE_DIR}/Code/data/solar/model_agss09.dat
            ${CMAKE_CURRENT_SOURCE_DIR}/Code/data/solar/Nuclear_Data.txt)
  set_tests_properties(single_collision PROPERTIES TIMEOUT 30)

  find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter)
  add_test(NAME project_contract
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Code/tests/test_project_contract.py)
  set_tests_properties(project_contract PROPERTIES TIMEOUT 10)
  add_test(NAME baseline_contract
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Code/tests/test_baseline.py)
  set_tests_properties(baseline_contract PROPERTIES TIMEOUT 10)
  add_test(NAME reference_physics_provenance
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/Code/tests/test_reference_physics_provenance.py)
  set_tests_properties(reference_physics_provenance PROPERTIES TIMEOUT 10)
  add_test(NAME t03_oracle_contract
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/Code/tests/test_t03_oracle_contract.py)
  set_tests_properties(t03_oracle_contract PROPERTIES TIMEOUT 10)

  add_test(NAME t03_validation_contract
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/Code/tests/test_t03_validation_runner.py
            $<TARGET_FILE:t03_validation_sampler>)
  set_tests_properties(t03_validation_contract PROPERTIES TIMEOUT 20)

  add_test(NAME grid_schema
    COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/Code/tests/test_grid_schema.py
            $<TARGET_FILE:t04_grid_schema_probe>)
  set_tests_properties(grid_schema PROPERTIES TIMEOUT 10)

  set_tests_properties(
    phase_space_grid
    phase_space_projector
    phase_space_cell_geometry
    collision_kernel
    scattering_physics
    solar_background
    scattering_rate
    target_sampling
    single_collision
    project_contract
    baseline_contract
    reference_physics_provenance
    t03_oracle_contract
    t03_validation_contract
    grid_schema
    PROPERTIES LABELS fast)
endif()
