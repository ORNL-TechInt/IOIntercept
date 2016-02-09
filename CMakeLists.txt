project(libiointercept)
cmake_minimum_required(VERSION 2.6)

set(CMAKE_C_COMPILER "cc")
set(CMAKE_C_FLAGS "-g")
add_subdirectory(src)

add_subdirectory(examples)
