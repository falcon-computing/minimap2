#ifndef MMSW_TASK_H
#define MMSW_TASK_H

#include <stdlib.h>
#include <stdexcept>
#include <time.h>

#include <glog/logging.h>

#include "blaze/Block.h" 
#include "blaze/Task.h" 

#include "blaze/xlnx_opencl/OpenCLEnv.h" 

#include <boost/thread/mutex.hpp>


class MnmpSWBuffer {
 public:
  MnmpSWBuffer(blaze::OpenCLEnv* env, int bank, size_t buffer_size) {
    //PLACE_TIMER;

    data = aligned_alloc(4096, buffer_size);
    memset(data, 0, buffer_size);

    // allocate input buffer
    static unsigned bankID[4] = {
      XCL_MEM_DDR_BANK0, XCL_MEM_DDR_BANK1,
      XCL_MEM_DDR_BANK2, XCL_MEM_DDR_BANK3};

    cl_mem_ext_ptr_t input_ext;

    input_ext.flags  = bankID[0];
    input_ext.obj    = (void*)data;
    input_ext.param  = 0;

    cl_int err = 0;
    cl_context context = env->getContext();
    buf = clCreateBuffer(context, 
        CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX | CL_MEM_USE_HOST_PTR,
        buffer_size, &input_ext, &err);

    if (!buf || err != CL_SUCCESS)
      throw blaze::invalidParam("failed to allocate CL buffer");

    DVLOG(1) << "create one mmsw buffer for bank: " << bank;
  }

  ~MnmpSWBuffer() {
    clReleaseMemObject(buf);
    free(data);  
    DVLOG(1) << "free one mmsw buffer";
  }
  
  cl_mem    buf;
  void*     data;
};

typedef boost::shared_ptr<MnmpSWBuffer> MnmpSWBuffer_ptr;

class MnmpSW : public blaze::Task {
 public:
  // extends the base class constructor
  // to indicate how many input blocks
  // are required
  MnmpSW();
  virtual ~MnmpSW();

  virtual uint64_t estimateClientTime() { return 60e9; }
  virtual uint64_t estimateTaskTime() { return 60e9; }

  void setupIndex(uint32_t*, uint32_t, uint64_t*, uint32_t*);
  virtual void compute();
  virtual void prepare();

 private:
  int conf_str2int(std::string key) {
    std::string conf;
    get_conf(key, conf);
    return std::stoi(conf);
  }

  // ksight::IntvTimer timer_;
  blaze::OpenCLEnv* env;

  static boost::mutex mtx_;
  static bool         index_init_;
  static uint32_t     num_ref_seq_;
  static cl_mem       opt_buffer_;
  static cl_mem       ref_buffer_;
  static cl_mem       mmi_buffer_;

  int                  cur_batch_size_;
  MnmpSWBuffer_ptr     input_;
  blaze::DataBlock_ptr output_;
};

// define the constructor and destructor for dlopen()
extern "C" blaze::Task* create() {
  return new MnmpSW();
}

extern "C" void destroy(blaze::Task* p) {
  delete p;
}

#endif