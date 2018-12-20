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
    int num_region): 
  client_(client),
  use_cpu_(false),
  input_(&alignInputs),
  num_region_(num_region)
{
  PLACE_TIMER;

  output_.resize(num_region);

  DLOG(INFO) << "fpga workload: num_region = " << num_region;
}

MnmpSWWorker:: ~MnmpSWWorker() {
  DLOG(INFO) << "Destroying worker";
}

// avx compute routine
void MnmpSWWorker::compute() {
  ksight::AutoTimer __timer("compute on cpu");

  if (use_cpu_) {
    mmsw_orig_compute(g_mnmpOpt, g_minimizer, *input_, output_);
  }
}

void MnmpSWWorker::getOutput(std::vector<align_output>& alignOutputs) {
  PLACE_TIMER;

  std::swap(alignOutputs, output_);
} 

// in this version, use a separate thread to run compute()
// for all reads/haps that is impossible for fpga to compute
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

    // each fpga client invokation takes a maximum data size
    {
      PLACE_TIMER1("serialize input");
      client_->setup(*input_, num_region_);
    }

    // start fpga run
    {
      PLACE_TIMER1("start()");
      client_->start();
    }

    // if cpu fallback is called, skip output copy
    {
      PLACE_TIMER1("copy output");
      
      // process output with data copy
      fpga_outputs *results = (fpga_outputs *)client_->getOutputPtr(0); 
      memcpy(output_.data(), results, sizeof(align_output)*num_region_);
    }
    // don't count cpu time
    ksight::ksight.add("compute on client", timer.stop());
    
    // t.join();
  }
}