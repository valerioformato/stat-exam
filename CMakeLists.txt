project(StatLab-Unfolding)
cmake_minimum_required(VERSION 3.13)

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

add_subdirectory(src)