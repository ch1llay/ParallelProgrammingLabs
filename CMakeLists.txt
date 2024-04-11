cmake_minimum_required(VERSION 3.25)
project(pp51)
project(pp52)

set(CMAKE_CXX_STANDARD 17)

# Add MPI Package to Project
find_package(MPI REQUIRED)

add_executable(pp51 one.cpp)
add_executable(pp52 two.cpp)

# Link MPI libraries
target_link_libraries(pp51 PRIVATE MPI::MPI_CXX)
target_link_libraries(pp52 PRIVATE MPI::MPI_CXX)