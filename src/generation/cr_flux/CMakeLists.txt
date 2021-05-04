set(CRSources main.cpp resolution_model.cpp spectrum.cpp)

add_executable(cr_generate_sample ${CRSources})
target_include_directories(cr_generate_sample PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(cr_generate_sample PRIVATE fmt::fmt-header-only ROOT::Core ROOT::RIO ROOT::Gpad ROOT::Tree ROOT::Hist)