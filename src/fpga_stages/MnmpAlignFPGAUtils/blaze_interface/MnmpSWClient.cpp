#include <glog/logging.h>
#include <vector>

#include "blaze/Client.h"
#include "MnmpSWClient.h"

using namespace blaze;

MnmpSWClient::MnmpSWClient(uint32_t* opt_data, uint32_t num_ref_seq, uint64_t* ref_data, uint32_t* mmi_data): 
  blaze::Client("MnmpSW", 6, 1)
{
  // skip input 0 for num_regions
  createInput(1, 1, sizeof(fpga_inputs), 1); // fpga_inputs

  createInput(2, 1, sizeof(mm_mapopt_t_union), 1); // options
  memcpy(getInputPtr(2), opt_data, sizeof(mm_mapopt_t_union));
  num_ref_seq_ = num_ref_seq;
  setInput(3, &num_ref_seq_, 1, 1, sizeof(uint32_t));
  createInput(4, 1, sizeof(uint64_t) * num_ref_seq_ * 2, 1); // references
  memcpy(getInputPtr(4), ref_data, sizeof(uint64_t)*num_ref_seq_*2);
  createInput(5, 1, sizeof(uint32_t*), 1); // minimizers
  memcpy(getInputPtr(5), &mmi_data, sizeof(uint32_t*));
  // memcpy(getInputPtr(5), mmi_data, MAX_REF_SIZE);
}


void MnmpSWClient::setup(
    std::vector<align_input>& inputs,
    int num_region) 
{
  PLACE_TIMER;
  num_region_ = num_region;

  uint64_t total_rl = 0;
  uint64_t total_hl = 0;

  setInput(0, &num_region_, 1, 1, sizeof(int));
  memcpy(getInputPtr(1), inputs.data(), sizeof(align_input)*num_region_);
  // memcpy(getInputPtr(1), inputs.data(), sizeof(fpga_inputs));
}

// load balance compute routine
void MnmpSWClient::compute() {
  ksight::AutoTimer __timer("compute on client cpu");
  ksight::ksight.add("num regions on cpu", num_region_);
}
