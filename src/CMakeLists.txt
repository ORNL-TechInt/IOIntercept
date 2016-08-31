include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/libemulatebbapi)
include_directories(/opt/ibm)
#Dynamic Target
set(CMAKE_C_COMPILER "gcc")
add_library(iointercept_dynamic SHARED libintercept.c)
SET_TARGET_PROPERTIES(iointercept_dynamic PROPERTIES COMPILE_DEFINITIONS "LIBINTERCEPT_PRELOAD")
#SET_TARGET_PROPERTIES(iointercept_dynamic PROPERTIES LINK_FLAGS "-L/opt/ibm/bb/lib -Wl,--rpath=/opt/ibm/bb/lib -lbbAPI")
#Static Target
add_library(iointercept libintercept.c)
#SET_TARGET_PROPERTIES(iointercept PROPERTIES LINK_FLAGS "-L/opt/ibm/bb/lib -Wl,--rpath=/opt/ibm/bb/lib -lbbAPI")
