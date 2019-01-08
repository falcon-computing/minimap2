set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wc++-compat -DHAVE_KALLOC")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wno-unused-function -O2")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -Wall")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB_RECURSE SRC_LIST "*.c")
file(GLOB_RECURSE KSW2_EXT_SRC_LIST "${CMAKE_CURRENT_SOURCE_DIR}/ksw2_ext*2_sse.c")

list(REMOVE_ITEM SRC_LIST "${CMAKE_CURRENT_SOURCE_DIR}/main.c"
                          "${CMAKE_CURRENT_SOURCE_DIR}/example.c"
                          "${CMAKE_CURRENT_SOURCE_DIR}/ksw2_ll_sse.c"
                          "${CMAKE_CURRENT_SOURCE_DIR}/ksw2_dispatch.c"
                          ${KSW2_EXT_SRC_LIST})


add_library(minimap-c-ksw2-ll STATIC "${CMAKE_CURRENT_SOURCE_DIR}/ksw2_ll_sse.c")
set_target_properties(minimap-c-ksw2-ll PROPERTIES
                      COMPILE_FLAGS "-msse2")

add_library(minimap-c-ksw2-extx-sse4 STATIC ${KSW2_EXT_SRC_LIST})
set_target_properties(minimap-c-ksw2-extx-sse4 PROPERTIES
                      COMPILE_FLAGS "-msse4.1 -DKSW_CPU_DISPATCH")

add_library(minimap-c-ksw2-extx-sse2 STATIC ${KSW2_EXT_SRC_LIST})
set_target_properties(minimap-c-ksw2-extx-sse2 PROPERTIES
                      COMPILE_FLAGS "-msse2 -mno-sse4.1 -DKSW_CPU_DISPATCH -DKSW_SSE2_ONLY")

add_library(minimap-c-ksw2-dispatch STATIC "${CMAKE_CURRENT_SOURCE_DIR}/ksw2_dispatch.c")
set_target_properties(minimap-c-ksw2-dispatch PROPERTIES
                      COMPILE_FLAGS "-msse4.1 -DKSW_CPU_DISPATCH")

add_library(minimap-c STATIC ${SRC_LIST})
target_link_libraries(minimap-c
  minimap-c-ksw2-ll
  minimap-c-ksw2-dispatch
  minimap-c-ksw2-extx-sse4
  minimap-c-ksw2-extx-sse2)
set_target_properties(minimap-c PROPERTIES OUTPUT_NAME minimap2)
