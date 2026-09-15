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

add_library(transport_reference_physics
  Code/src/physics/ScatteringPhysics.cpp
  Code/src/physics/SolarBackground.cpp)
target_include_directories(transport_reference_physics PUBLIC Code/include)
target_compile_features(transport_reference_physics PUBLIC cxx_std_14)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(transport_reference_physics PRIVATE
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
endif()
