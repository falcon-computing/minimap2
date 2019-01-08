add_custom_target(Xilinx)

if(DEFINED ENV{XILINX_SDX})
  set(Xilinx_DIR "$ENV{XILINX_SDX}")
  set(Xilinx_INCLUDE_DIRS "${Xilinx_DIR}/runtime/include/1_2")
  set(Xilinx_LIBRARY_DIRS "${Xilinx_DIR}/runtime/lib/x86_64")
  set(Xilinx_LIBRARIES "xilinxopencl")
  set(Xilinx_FOUND TRUE)

  message(STATUS "Found Xilinx OpenCL Environment")
endif() 