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

include(CTest)
if(BUILD_TESTING)
  add_executable(test_phase_space_grid Code/tests/test_phase_space_grid.cpp)
  target_link_libraries(test_phase_space_grid PRIVATE transport_grid)
  add_test(NAME phase_space_grid COMMAND test_phase_space_grid)
  set_tests_properties(phase_space_grid PROPERTIES TIMEOUT 10)

  find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter)
  add_test(NAME project_contract
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Code/tests/test_project_contract.py)
  set_tests_properties(project_contract PROPERTIES TIMEOUT 10)
  add_test(NAME baseline_contract
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Code/tests/test_baseline.py)
  set_tests_properties(baseline_contract PROPERTIES TIMEOUT 10)
endif()
