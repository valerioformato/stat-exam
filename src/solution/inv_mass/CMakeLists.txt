set(IMSources main.cpp)

add_executable(im_solution ${IMSources})
target_include_directories(im_solution PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(im_solution PRIVATE fmt::fmt-header-only ROOT::Core ROOT::RIO ROOT::Gpad ROOT::Tree ROOT::Hist)