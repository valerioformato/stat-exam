project(StatLab-Unfolding)
cmake_minimum_required(VERSION 3.13)

# Check for C++ standard and set it to C++14 if not already set
set(CMAKE_CXX_STANDARD_REQUIRED 14)
set(CMAKE_CXX_STANDARD 14)

# handle dependencies
include(FetchContent)

FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 6.1.2)
FetchContent_GetProperties(fmt)
if(NOT fmt_POPULATED)
  set(FMT_INSTALL ON)
  FetchContent_Populate(fmt)
  add_subdirectory(${fmt_SOURCE_DIR} ${fmt_BINARY_DIR})
endif()

find_package(ROOT CONFIG REQUIRED)

if(DEFINED ASAN)
  add_compile_options(-fsanitize=address -fsanitize=leak)
  add_link_options(-fsanitize=address -fsanitize=leak)
endif()

add_subdirectory(src)