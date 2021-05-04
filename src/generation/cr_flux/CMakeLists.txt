set(CRSources main.cpp resolution_model.cpp spectrum.cpp)

add_executable(generate_sample ${CRSources})
target_include_directories(generate_sample PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(generate_sample PRIVATE fmt::fmt-header-only ROOT::Core ROOT::RIO ROOT::Gpad ROOT::Tree ROOT::Hist)