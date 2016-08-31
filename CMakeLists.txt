project(libiointercept)
cmake_minimum_required(VERSION 2.6)


option(DEBUG "Option description" OFF)



set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
set(CMAKE_C_FLAGS "-g -O2")

if(DEBUG)
  add_definitions(-DDEBUG)
endif()

add_subdirectory(src)

add_subdirectory(libemulatebbapi)

add_subdirectory(examples)

