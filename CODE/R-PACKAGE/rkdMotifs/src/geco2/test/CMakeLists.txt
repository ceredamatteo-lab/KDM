cmake_minimum_required(VERSION 3.0)

project(test_geco2)

find_package(geco2 REQUIRED)
find_package(OpenMP REQUIRED)

add_executable(test_bed test_beds.cpp)

target_link_libraries(test_bed PRIVATE geco2::geco2 OpenMP::OpenMP_CXX)



add_executable(test_counts test_counts.cpp)

target_link_libraries(test_counts PRIVATE geco2::geco2 OpenMP::OpenMP_CXX)
