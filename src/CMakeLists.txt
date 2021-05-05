configure_file(localization.h.in ${CMAKE_CURRENT_BINARY_DIR}/localization.h @ONLY)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(generation)
add_subdirectory(solution)

