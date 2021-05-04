set(CRSources main.cpp resolution_model.cpp)

add_executable(im_generate_sample ${CRSources})
target_include_directories(im_generate_sample PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(im_generate_sample PRIVATE fmt::fmt-header-only ROOT::Core ROOT::RIO ROOT::Gpad ROOT::Tree ROOT::Hist ROOT::Physics)