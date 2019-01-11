#ifndef MMSWWORKER_H
#define MMSWWORKER_H

#include <vector>
#include <stdexcept>

#include "blaze/Client.h"
#include "MnmpSWClient.h"
#include "../MnmpAlignFPGAUtils/host_interface/host_types.h"

class MnmpSWWorker {
 public:
  MnmpSWWorker(MnmpSWClient* client,
    std::vector<align_input>& alignInputs,
    std::vector<align_output>& alignOutputs,
    int num_region);

  ~MnmpSWWorker();

  // perform all computation
  void run();

  void getOutput(std::vector<align_output>& alignOutputs);

  void compute();     // fallback compute() when client has problem

 private:
  void distribute() {;} // TODO: distribute workloads using wait time

  MnmpSWClient* client_;

  bool use_cpu_;      // if true, use cpu compute everything
  int num_region_;

  std::vector<align_input>  *input_;
#ifdef DUPLICATE_OUTPUT
  std::vector<align_output>  output_;
#else
  std::vector<align_output> *output_;
#endif

  uint64_t cpu_end_ts_;  // finish timestamp of cpu compute
  uint64_t fpga_end_ts_; // finish timestamp of fpga compute
};


#endif
