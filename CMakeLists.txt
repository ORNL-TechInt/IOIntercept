project(libiointercept)
cmake_minimum_required(VERSION 2.6)


option(DEBUG "Use debugging symbols" OFF)
option(TITAN "Build for titan emulation" OFF)


set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)

if(DEBUG)
  add_definitions(-DDEBUG)
  set(CMAKE_C_FLAGS "-g -Werror")
else()
  set(CMAKE_C_FLAGS "-g -O3 -Werror")
endif()

add_subdirectory(src)

if (TITAN)
add_definitions(-DTITAN)
endif()

add_subdirectory(examples)

file(COPY test DESTINATION .)
