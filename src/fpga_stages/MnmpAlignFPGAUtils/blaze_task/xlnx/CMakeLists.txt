include_directories(${Xilinx_INCLUDE_DIRS})

link_directories(${Xilinx_LIBRARY_DIRS})

set(proj_name MnmpSWTask_xlnx)
add_library(${proj_name} SHARED MnmpSWTask.cpp)
add_dependencies(${proj_name} Blaze)
add_dependencies(${proj_name} Google)
add_dependencies(${proj_name} KSight)

target_link_libraries(${proj_name}
  blaze
  blaze_message
  ${Xilinx_LIBRARIES})

set_target_properties(${proj_name}
  PROPERTIES OUTPUT_NAME MnmpSWTask)

install(TARGETS MnmpSWTask_xlnx
    LIBRARY DESTINATION mmsw_task/xlnx)