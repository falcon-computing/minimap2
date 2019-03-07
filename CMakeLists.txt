include(ExternalProject)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

cmake_minimum_required(VERSION 2.8.12)

project(minimap-flow)

# set up dependency
find_package(Boost 1.53.0 REQUIRED
             COMPONENTS system thread chrono iostreams filesystem regex)
find_package(ZLIB     REQUIRED)
find_package(Threads  REQUIRED)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(Google   REQUIRED)
find_package(FalconLM REQUIRED)
find_package(HTS      REQUIRED)
find_package(Xilinx)

# build option
option( USE_FPGA   "Build with FPGA"                        OFF )
option( USE_BLAZE  "Build with Blaze"                       ON )
option( USE_NATIVE "Build with native FPGA routines"        OFF)
option( USE_CKERNEL   "Build with original C code for kernels" OFF)
# priority: USE_BLAZE < USE_NATIVE < USE_CKERNEL

if (XilinxOCL_FOUND AND USE_FPGA)
 set(FPGA_TYPE "xlnx" )
 set(FPGA_INCLUDE_DIRS "${Xilinx_INCLUDE_DIRS}" )
 set(FPGA_LIBRARY_DIRS "${Xilinx_LIBRARY_DIRS}")
 set(FPGA_LIBRARIES "${Xilinx_LIBRARIES}")
#elseif ( DEFINED ENV{ALTERAOCLSDKROOT} AND USE_FPGA)
#  set( USE_FPGA On )
#  set( FPGA_TYPE "intel" )
#  find_package( IntelAltera REQUIRED )
#  set( FPGA_INCLUDE_DIRS "${IntelAltera_INCLUDE_DIRS}" )
#  set( FPGA_LIBRARY_DIRS "${IntelAltera_LIBRARY_DIRS}")
#  set( FPGA_LIBRARIES "${IntelAltera_LIBRARIES}")
else()
  message( STATUS "No OpenCL module loaded. Build without using FPGA." )
  set( USE_FPGA Off )
endif()

if ( USE_FPGA )
  # find_package( Protobuf       REQUIRED )
  find_package( Blaze          REQUIRED )
  find_package( KSight         REQUIRED )
endif()


# set up debug/release mode
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Debug")
  set(CMAKE_BUILD_TYPE "Debug")
endif()
message( STATUS "Build type: ${CMAKE_BUILD_TYPE}" )

add_subdirectory(kflow)
add_subdirectory(minimap2)

# set up for compilation
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -fPIC -DUSE_HTSLIB -msse4.1")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -rdynamic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -rdynamic -Wall -DUSELICENSE")
if ( USE_MPI )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_LINK_FLAGS}")
endif()


# get code version
execute_process( COMMAND git describe --tags
  COMMAND  sed "s/\\(.*\\)-.*/\\1/"
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE )
set(BIN_VERSION "${GIT_VERSION}")


if ( USE_FPGA )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBUILD_FPGA" )
  if ( USE_CKERNEL )
    message( STATUS "Build with original C kernel for FPGA" )
    set(USE_NATIVE OFF)
    set(USE_BLAZE OFF)
  elseif ( USE_NATIVE )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_NATIVE" )
    message( STATUS "Build with OpenCL routine for FPGA" )
    set(USE_BLAZE OFF)
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLOCAL_BLAZE" )
    message( STATUS "Build with Blaze library for FPGA" )
  endif()

  set(BIN_VERSION "${BIN_VERSION}-${FPGA_TYPE}")
  if ( ${FPGA_TYPE} MATCHES "xlnx" )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DXILINX_FPGA" )
  elseif ( ${FPGA_TYPE} MATCHES "intel" )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DINTEL_FPGA" )
  endif()
endif()

if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  set(BIN_VERSION "${BIN_VERSION}-dev")
endif()
if (DEFINED DEPLOYMENT_DST)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEPLOY_${DEPLOYMENT_DST}")
  set(BIN_VERSION "${BIN_VERSION}-${DEPLOYMENT_DST}")
endif()

message(STATUS "version: ${BIN_VERSION}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERSION=\\\"${BIN_VERSION}\\\"")

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/src/cpu_stages)

include_directories( ./kflow/include )
include_directories( ./minimap2 )
link_directories( ./minimap2 )

include_directories(
  ${Google_INCLUDE_DIRS}
  ${FalconLM_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${HTS_INCLUDE_DIRS})

link_directories(
  ${Google_LIBRARY_DIRS}
  ${FalconLM_LIBRARY_DIRS}
  ${Boost_LIBRARY_DIRS}
  ${HTS_LIBRARY_DIRS})

file(GLOB MM_SRC_LIST "src/*.cpp" "src/cpu_stages/*.cpp")
list(REMOVE_ITEM MM_SRC_LIST "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")

if ( USE_FPGA )
  include(${CMAKE_CURRENT_SOURCE_DIR}/src/fpga_stages/CMakeLists.txt)
endif()


add_executable(minimap-flow ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp ${MM_SRC_LIST})
if ( USE_FPGA )
  add_dependencies(minimap-flow Blaze KSight)
endif()
add_dependencies(minimap-flow Google FalconLM HTS)

set (CMAKE_STATIC_LINKER_FLAGS "--as-needed")

if ( USE_FPGA )
  target_link_libraries(minimap-flow
    minimap-c
    kflow
    numa
    ${Boost_LIBRARIES} 
    ${FalconLM_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${ZLIB_LIBRARIES}
    ${CMAKE_DL_LIBS}
    ${HTS_LIBRARIES}
    ${Blaze_LIBRARIES}
    ${KSight_LIBRARIES}
    ${FPGA_LIBRARIES}
    ${Google_LIBRARIES})
else()
  target_link_libraries(minimap-flow
    minimap-c
    kflow 
    numa
    ${Boost_LIBRARIES} 
    ${FalconLM_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${ZLIB_LIBRARIES}
    ${CMAKE_DL_LIBS}
    ${HTS_LIBRARIES}
    ${Google_LIBRARIES})
endif()

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
   message(STATUS "Install prefix is set to ~/.falcon-genome/minimap/${BIN_VERSION}")
   set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.falcon-genome/minimap/${BIN_VERSION}")
endif()
install(TARGETS minimap-flow RUNTIME DESTINATION .)

enable_testing()
add_subdirectory(test)
