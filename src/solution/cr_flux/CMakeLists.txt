set(CRSources main.cpp)

add_executable(cr_solution ${CRSources})
target_include_directories(cr_solution PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(cr_solution PRIVATE fmt::fmt-header-only ROOT::Core ROOT::RIO ROOT::Gpad ROOT::Tree ROOT::Hist)