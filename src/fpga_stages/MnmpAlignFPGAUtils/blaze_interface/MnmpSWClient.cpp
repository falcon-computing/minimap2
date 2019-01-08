#include <glog/logging.h>
#include <vector>

#include "blaze/Client.h"
#include "MnmpSWClient.h"

using namespace blaze;

MnmpSWClient::MnmpSWClient(uint32_t* opt_data, uint32_t num_ref_seq, uint64_t* ref_data, uint32_t* mmi_data): 
  blaze::Client("MnmpSW", NUM_INPUT_ARGS, 1)
{
  opt_data_    = opt_data;
  num_ref_seq_ = num_ref_seq;
  ref_data_    = ref_data;
  mmi_data_    = mmi_data;

  setInput(2, &opt_data_,    1, 1, sizeof(uint32_t*));
  setInput(3, &num_ref_seq_, 1, 1, sizeof(uint32_t));
  setInput(4, &ref_data_,    1, 1, sizeof(uint64_t*));
  setInput(5, &mmi_data,     1, 1, sizeof(uint32_t*));
}


void MnmpSWClient::setup(
    std::vector<align_input>& inputs,
#ifndef DUPLICATE_OUTPUT
    std::vector<align_output>& outputs,
#endif
    int num_region) 
{
  PLACE_TIMER;
  num_region_ = num_region;
  inputs_data_ = (fpga_inputs*)inputs.data();

  setInput(0, &num_region_,  1, 1, sizeof(int));
  setInput(1, &inputs_data_, 1, 1, sizeof(fpga_inputs*));
  // memcpy(getInputPtr(1), inputs.data(), sizeof(align_input)*num_region_);

#ifndef DUPLICATE_OUTPUT
  outputs_data_ = (fpga_outputs*)outputs.data();
  setInput(6, &outputs_data_, 1, 1, sizeof(fpga_outputs*));
#endif
}

// load balance compute routine
void MnmpSWClient::compute() {
  ksight::AutoTimer __timer("compute on client cpu");
  ksight::ksight.add("num regions on cpu", num_region_);
}
