cmake_minimum_required(VERSION 3.22.2)

project(main)

find_package(MPI REQUIRED)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(main MPI::MPI_CXX)