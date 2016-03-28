add_library(emulate_bbapi_dynamic SHARED bbapi.c bbapi.h)
add_library(emulate_bbapi  bbapi.c bbpublic.c bbapi.h)
set (LIBEMULATEBBAPI_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

