include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/libemulatebbapi)
#Dynamic Target
add_library(iointercept_dynamic SHARED libintercept.c)
SET_TARGET_PROPERTIES(iointercept_dynamic PROPERTIES COMPILE_DEFINITIONS "LIBINTERCEPT_PRELOAD")

#Static Target
add_library(iointercept libintercept.c)
