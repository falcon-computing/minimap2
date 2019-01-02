#include <glog/logging.h>
#include <unordered_set>
#include <vector>

#include "blaze/Client.h"
#include "ksight/tools.h"
#include "MnmpSWWorker.h"
#include "../MnmpAlignFPGAUtils/host_interface/mmsw_orig.h"
#include "../MnmpAlignFPGAUtils/host_interface/mmsw_fpga.h"

#include "MnmpGlobal.h"

using namespace blaze;

MnmpSWWorker::MnmpSWWorker(
    MnmpSWClient* client,
    std::vector<align_input>& alignInputs,
    std::vector<align_output>& alignOutputs,
    int num_region): 
  client_(client),
  use_cpu_(false),
  input_(&alignInputs),
#ifndef DUPLICATE_OUTPUT
  output_(&alignOutputs),
#endif
  num_region_(num_region)
{
  PLACE_TIMER;

#ifdef DUPLICATE_OUTPUT
  output_.resize(num_region);
#else
  output_->resize(num_region);
#endif

  DLOG(INFO) << "fpga workload: num_region = " << num_region;
}

MnmpSWWorker:: ~MnmpSWWorker() {
  DLOG(INFO) << "Destroying worker";
}

// cpu compute routine
void MnmpSWWorker::compute() {
  ksight::AutoTimer __timer("compute on cpu");

  if (use_cpu_) {
#ifdef DUPLICATE_OUTPUT
    mmsw_orig_compute(g_mnmpOpt, g_minimizer, *input_, output_);
#else
    mmsw_orig_compute(g_mnmpOpt, g_minimizer, *input_, *output_);
#endif
  }
}

void MnmpSWWorker::getOutput(std::vector<align_output>& alignOutputs) {
  PLACE_TIMER;

#ifdef DUPLICATE_OUTPUT
  // swap the buffers between worker and compute_stage
  std::swap(alignOutputs, output_);
#endif
} 

// currently all workload are calculated on fpga
void MnmpSWWorker::run() {
  DLOG(INFO) << "Computing mode: " << (use_cpu_ ? "CPU" : "CFX");
  if (use_cpu_) {
    compute();
  }
  else {
    //PLACE_TIMER1("compute on fpga");
    ksight::Timer timer;
    timer.start();

    // start a thread to run cpu
    // boost::thread t(boost::bind(&PairHMMWorker::compute, this));

    // set input
    {
      PLACE_TIMER1("serialize input");
#ifdef DUPLICATE_OUTPUT
      client_->setup(*input_, num_region_);
#else
      client_->setup(*input_, *output_, num_region_);
#endif
    }

    // start fpga run
    {
      PLACE_TIMER1("start()");
      client_->start();
    }

    // if cpu fallback is called, skip output copy
    {
      PLACE_TIMER1("copy output");

#ifdef DUPLICATE_OUTPUT
      // copy the data from MMSWTask output buffer
      fpga_outputs *results = (fpga_outputs *)client_->getOutputPtr(0); 
      memcpy(output_.data(), results, sizeof(align_output)*num_region_);
#else
      // get a pseudo-output token
      int ret = *((int *)client_->getOutputPtr(0));
      // assert(ret == 1);
#endif
    }
    // don't count cpu time
    ksight::ksight.add("compute on client", timer.stop());
    
    // t.join();
  }
}