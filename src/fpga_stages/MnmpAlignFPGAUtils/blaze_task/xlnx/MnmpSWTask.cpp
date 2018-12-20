#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <time.h>

#include <glog/logging.h>

// #include "ksight/tools.h"
#include "MnmpSWTask.h"

#include "blaze/xlnx_opencl/OpenCLEnv.h" 

#include "../MnmpAlignFPGAUtils/host_interface/host_types.h"
#include "../MnmpAlignFPGAUtils/host_interface/mmsw_fpga.h"
// #include "mmsw.h"

boost::mutex MnmpSW::mtx_;

bool      MnmpSW::index_init_ = false;
uint32_t  MnmpSW::num_ref_seq_ = 0;
cl_mem    MnmpSW::opt_buffer_ = NULL;
cl_mem    MnmpSW::ref_buffer_ = NULL;
cl_mem    MnmpSW::mmi_buffer_ = NULL;


MnmpSW::MnmpSW(): blaze::Task(6)//, timer_("MnmpSW Task")
{
  ;
}

MnmpSW::~MnmpSW() {
  if (env) {
    // push back the scratch blocks to TaskEnv
    env->putScratch("input",  input_);
    env->putScratch("output", output_);
  }
}

void MnmpSW::setupIndex(uint32_t* opt_data, uint32_t num_ref_seq, uint64_t* ref_data, uint32_t* mmi_data) {
  mtx_.lock();
  if (index_init_) {
    mtx_.unlock();
    return;
  }

  int opt_bank    = 0;
  int mmi_bank    = 0;
  int ref_bank    = 0;

  unsigned bankID[4] = {XCL_MEM_DDR_BANK0, XCL_MEM_DDR_BANK1, \
        XCL_MEM_DDR_BANK2, XCL_MEM_DDR_BANK3};
  cl_mem_ext_ptr_t bank_ext[4];
  bank_ext[0].flags = bankID[0];
  bank_ext[0].param = 0;
  bank_ext[0].obj = nullptr;
  cl_int err = 0;
  cl_int cl_status;
  cl_context context = env->getContext();
  cl_command_queue command = env->getCmdQueue();

  {
    //PLACE_TIMER1("init option buffer");
    opt_buffer_ = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
                      sizeof(mm_mapopt_t_union), &bank_ext[0], &err);
    if (!opt_buffer_)
        printf("failed to create opt_buffer, err code is %d\n", err);
    else 
        printf("create opt_buffer\n");
    
    cl_status = clEnqueueWriteBuffer(command, opt_buffer_, CL_TRUE,\
                  0, sizeof(mm_mapopt_t_union), opt_data, 0, NULL, NULL);
    if (cl_status != CL_SUCCESS) {
        printf("failed to write opt buffer, err code is %d\n", err);
    }
  }

  {
    //PLACE_TIMER1("init reference buffer");
    ref_buffer_ = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
                      sizeof(uint64_t) * num_ref_seq * 2, &bank_ext[0], &err);
    if (!ref_buffer_)
        printf("failed to create ref_buffer, err code is %d\n", err);
    else 
        printf("create ref_buffer\n");
    
    cl_status = clEnqueueWriteBuffer(command, ref_buffer_, CL_TRUE,\
                    0, sizeof(uint64_t) * num_ref_seq * 2, ref_data , 0, NULL, NULL);
    if (cl_status != CL_SUCCESS) {
        printf("failed to write ref buffer, err code is %d\n", err);
    }
  }

  {
    //PLACE_TIMER1("init minimizer buffer");
    mmi_buffer_ = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
                      MAX_REF_SIZE, &bank_ext[0], &err);
    if (!mmi_buffer_)
        printf("failed to create mmi_buffer, err code is %d\n", err);
    else 
        printf("create mmi_buffer\n");

    cl_status = clEnqueueWriteBuffer(command, mmi_buffer_, CL_TRUE,\
                    0, MAX_REF_SIZE, mmi_data , 0, NULL, NULL);
    if (cl_status != CL_SUCCESS) {
        printf("failed to write mmi buffer, err code is %d\n", err);
    }
  }


  index_init_ = true;
  mtx_.unlock();
}

void MnmpSW::prepare() {
  PLACE_TIMER;
  
  env = (blaze::OpenCLEnv*)getEnv();

  cur_batch_size_ = *((int*)getInput(0));
  fpga_inputs *in_data = (fpga_inputs *)getInput(1);

  uint32_t *opt_data = (uint32_t *)getInput(2);
  num_ref_seq_ = *((uint32_t*)getInput(3));
  uint64_t *ref_data = (uint64_t *)getInput(4);
  uint32_t *mmi_data = *((uint32_t **)getInput(5));

  // getting configurations and prepare fpga input
  int input_bank  = 0;  
  int output_bank = 0;

  // init index if necessary
  if (!index_init_) {
    setupIndex(opt_data, num_ref_seq_, ref_data, mmi_data);
  }

  // get input block from scratch, create if it does not exists
  {
    //PLACE_TIMER1("create input buffer");
    if (!env->getScratch("input", input_)) {
      MnmpSWBuffer_ptr input(new MnmpSWBuffer(env, input_bank, sizeof(fpga_inputs)));
      input_ = input;
      DLOG(INFO) << "Creating a new MnmpSWBuffer for input";
    }
    else {
      DLOG(INFO) << "Getting a MnmpSWBuffer for input from Scratch";
    }
  }

  // allocate output buffer
  blaze::ConfigTable_ptr output_conf(new blaze::ConfigTable());
  output_conf->write_conf("bankID", output_bank);

  {
    //PLACE_TIMER1("create output buffer");
    if (!env->getScratch("output", output_)) {
      blaze::DataBlock_ptr b = env->create_block(1, 
          1, sizeof(fpga_outputs), 
          0, blaze::DataBlock::OWNED, output_conf);

      output_ = b;
    }
    setOutput(0, output_);
  }

  // write input buffer
  cl_int err = 0;
  cl_command_queue command = env->getCmdQueue();

  if (!command) {
    DLOG(ERROR) << "failed to get OpenCLEnv";
    throw blaze::invalidParam(__func__);
  }
  {
    //PLACE_TIMER1("write_buffer");
    cl_event event;
    cl_int err = clEnqueueWriteBuffer(command, input_->buf, CL_TRUE,\
                    0, sizeof(align_input) * cur_batch_size_, (uint32_t*)in_data , 0, NULL, &event);
    if (err != CL_SUCCESS) {
      throw std::runtime_error("failed to migrate input buffer");
    }
    clWaitForEvents(1, &event);
    clReleaseEvent(event);
  }
}

void MnmpSW::compute() {
  std::string kernel_name;
  get_conf("kernel_name", kernel_name);
  // ksight::IntvTimer i_timer("mmsw compute");

  //PLACE_TIMER1(kernel_name);
  uint64_t start_ts = blaze::getUs();

  cl_kernel        kernel  = env->getKernel();
  cl_command_queue command = env->getCmdQueue();

  int err = 0;
  cl_event event;
  err |= clSetKernelArg(kernel, 0, sizeof(cl_mem), &opt_buffer_);
  err |= clSetKernelArg(kernel, 1, sizeof(uint32_t), &num_ref_seq_);
  err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &ref_buffer_);
  err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &mmi_buffer_);
  err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &input_->buf);
  err |= clSetKernelArg(kernel, 5, sizeof(cl_mem), (cl_mem*)output_->getData());
  err |= clSetKernelArg(kernel, 6, sizeof(int), &cur_batch_size_);

  // uint64_t w_start = blaze::getUs();

  { //PLACE_TIMER1(kernel_name + " kernel");
  cl_event event;
  cl_int err = clEnqueueTask(command, kernel, 0, NULL, &event);
  if (err != CL_SUCCESS) {
    throw std::runtime_error("failed to enqueue kernel");
  }
  clWaitForEvents(1, &event);

// #ifndef NO_PROFILE
//   cl_ulong k_start, k_end;
//   clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(k_start), &k_start, NULL);
//   clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(k_end), &k_end, NULL);

//   ksight::ksight.add(kernel_name + " cells", (uint64_t)num_cell);
//   ksight::Interval<uint64_t> intv(start_ts, blaze::getUs());
//   ksight::ksight.add("mmsw kernel", intv);
// #endif

  clReleaseEvent(event);
  }
}