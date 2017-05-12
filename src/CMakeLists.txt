include_directories(${CMAKE_SOURCE_DIR}/include)
#include_directories(/opt/ibm)
include_directories(/usr/local/)

#Dynamic Target
add_library(spectral SHARED spectral.c util.c emapi.c empublic.c wrappers.c)
SET_TARGET_PROPERTIES(spectral PROPERTIES COMPILE_DEFINITIONS "SPECTRAL_PRELOAD")
SET_TARGET_PROPERTIES(spectral PROPERTIES LINK_FLAGS "-Wl,--rpath /usr/local/bb/lib -L/usr/local/bb/lib -lbbAPI")
#Static Target
add_library(spectralstatic spectral.c util.c emapi.c empublic.c wrappers.c)
set_target_properties(spectralstatic PROPERTIES OUTPUT_NAME spectral)

install(TARGETS spectral DESTINATION lib)
install(TARGETS spectralstatic DESTINATION lib)
