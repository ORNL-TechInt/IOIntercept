find_package(MPI)
set(CMAKE_C_COMPILER "mpicc")

if(DEBUG)
  add_definitions(-DDEBUG)
  set(CMAKE_C_FLAGS "-g -Werror -std=c99")
else()
  set(CMAKE_C_FLAGS "-g -O2 -Werror -std=c99")
endif()

if (TITAN)
add_definitions(-DTITAN)
endif()



add_executable(test_dynamic test.c)
add_executable(test_static test.c)
add_executable(test_mpi_static test_mpi.c)

if (TITAN)
include_directories(${CMAKE_SOURCE_DIR}/libemulatebbapi)
endif()

target_link_libraries(test_static spectralstatic)
target_link_libraries(test_mpi_static spectralstatic)

if (TITAN)
target_link_libraries(test_static emulate_bbapi)
target_link_libraries(test_mpi_static emulate_bbapi)
target_link_libraries(test_dynamic emulate_bbapi_dynamic)

set_target_properties(test_dynamic PROPERTIES LINK_FLAGS "-dynamic ")
set_target_properties(test_mpi_static PROPERTIES LINK_FLAGS "-Wl,--wrap,close -Wl,--wrap,fclose ")
set_target_properties(test_static PROPERTIES LINK_FLAGS "-Wl,--wrap,close -Wl,--wrap,fclose ")
else()
set_target_properties(test_static PROPERTIES LINK_FLAGS "-Wl,--wrap,close -Wl,--wrap,fclose -Wl,--rpath /opt/ibm/bb/lib -L/opt/ibm/bb/lib -lbbAPI ")
set_target_properties(test_dynamic PROPERTIES LINK_FLAGS "-dynamic -ldl -Wl,--rpath /opt/ibm/bb/lib -L/opt/ibm/bb/lib -lbbAPI ")
set_target_properties(test_mpi_static PROPERTIES LINK_FLAGS "-Wl,--wrap,close -Wl,--wrap,fclose -Wl,--rpath /opt/ibm/bb/lib -L/opt/ibm/bb/lib -lbbAPI ")
endif()

file(COPY gtc DESTINATION .)
file(COPY setupenv DESTINATION .)
