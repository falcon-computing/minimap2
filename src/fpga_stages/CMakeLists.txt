include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/fpga_stages
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/fpga_stages/MnmpAlignFPGAUtils
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/fpga_stages/MnmpAlignFPGAUtils/blaze_interface
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/fpga_stages/MnmpAlignFPGAUtils/host_interface)

include_directories(${FPGA_INCLUDE_DIRS} ${Blaze_INCLUDE_DIRS} ${KSight_INCLUDE_DIRS})
link_directories(${FPGA_LIBRARY_DIRS} ${Blaze_LIBRARY_DIRS} ${KSight_LIBRARY_DIRS})

file(GLOB MM_FPGA_SRC_LIST "src/fpga_stages/*.cpp"
       "src/fpga_stages/MnmpAlignFPGAUtils/*.cpp"
       "src/fpga_stages/MnmpAlignFPGAUtils/blaze_interface/*.cpp"
       "src/fpga_stages/MnmpAlignFPGAUtils/host_interface/*.cpp")
list(APPEND MM_SRC_LIST ${MM_FPGA_SRC_LIST})

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src/fpga_stages/MnmpAlignFPGAUtils/blaze_task/${FPGA_TYPE})