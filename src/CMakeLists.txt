include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

if (TITAN)
include_directories(${CMAKE_SOURCE_DIR}/libemulatebbapi)
add_definitions(-DTITAN)
else()
include_directories(/opt/ibm)
endif()

#Dynamic Target
add_library(iointercept SHARED libintercept.c)
SET_TARGET_PROPERTIES(iointercept PROPERTIES COMPILE_DEFINITIONS "LIBINTERCEPT_PRELOAD")
SET_TARGET_PROPERTIES(iointercept PROPERTIES LINK_FLAGS "-Wl,--rpath /opt/ibm/bb/lib -L/opt/ibm/bb/lib -lbbAPI")
#Static Target
add_library(iointerceptstatic libintercept.c)
set_target_properties(iointerceptstatic PROPERTIES OUTPUT_NAME iointercept)

install(TARGETS iointercept DESTINATION lib)
install(TARGETS iointerceptstatic DESTINATION lib)
