cmake_minimum_required(VERSION 3.12)
project(ray)

include_directories(${CMAKE_SOURCE_DIR}/src)

add_subdirectory(src)

add_executable(ray main.c)
