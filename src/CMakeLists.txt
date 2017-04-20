if (TITAN)
include_directories(${CMAKE_SOURCE_DIR}/libemulatebbapi)
add_definitions(-DTITAN)
else()
include_directories(/opt/ibm)
endif()

#Dynamic Target
add_library(spectral SHARED libintercept.c)
SET_TARGET_PROPERTIES(spectral PROPERTIES COMPILE_DEFINITIONS "LIBINTERCEPT_PRELOAD")
SET_TARGET_PROPERTIES(spectral PROPERTIES LINK_FLAGS "-Wl,--rpath /opt/ibm/bb/lib -L/opt/ibm/bb/lib -lbbAPI")
#Static Target
add_library(spectralstatic libintercept.c)
set_target_properties(spectralstatic PROPERTIES OUTPUT_NAME spectral)

install(TARGETS spectral DESTINATION lib)
install(TARGETS spectralstatic DESTINATION lib)
