ExternalProject_Add(glog-download
    PREFIX "glog"
    URL ${DEPS}/glog-falcon.tar.gz
    SOURCE_DIR "${CMAKE_BINARY_DIR}/glog/install"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

ExternalProject_Add(gflags-download
    PREFIX "gflags"
    URL ${DEPS}/gflags.tar.gz
    SOURCE_DIR "${CMAKE_BINARY_DIR}/gflags/install"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

ExternalProject_Add(googletest-download
    PREFIX "googletest"
    URL ${DEPS}/googletest.tar.gz
    SOURCE_DIR "${CMAKE_BINARY_DIR}/googletest/install"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

ExternalProject_Add(protobuf-download
    PREFIX "protobuf"
    URL ${DEPS}/protobuf-2.5.0.tar.gz
    SOURCE_DIR "${CMAKE_BINARY_DIR}/protobuf/install"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")

add_custom_target(Google)
add_dependencies(Google
    glog-download
    gflags-download
    googletest-download
    protobuf-download)

set(Google_INCLUDE_DIRS
    ${CMAKE_BINARY_DIR}/glog/install/include
    ${CMAKE_BINARY_DIR}/gflags/install/include
    ${CMAKE_BINARY_DIR}/googletest/install/include
    ${CMAKE_BINARY_DIR}/gtest/install/include
    ${CMAKE_BINARY_DIR}/protobuf/install/include)

set(Google_LIBRARY_DIRS
    ${CMAKE_BINARY_DIR}/glog/install/lib
    ${CMAKE_BINARY_DIR}/gflags/install/lib
    ${CMAKE_BINARY_DIR}/googletest/install/lib
    ${CMAKE_BINARY_DIR}/protobuf/install/lib)

set(Google_LIBRARIES
    glog
    gflags
    protobuf
    gtest)
