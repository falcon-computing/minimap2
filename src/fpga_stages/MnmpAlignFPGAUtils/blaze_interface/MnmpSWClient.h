#ifndef MMSWCLIENT_H
#define MMSWCLIENT_H
#include <stdexcept>

#include "blaze/Client.h"
#include "ksight/tools.h"
#include "../MnmpAlignFPGAUtils/host_interface/host_types.h"
#include "../MnmpAlignFPGAUtils/host_interface/mmsw_fpga.h"

#ifdef DUPLICATE_OUTPUT
#define NUM_INPUT_ARGS (6)
#else
#define NUM_INPUT_ARGS (7)
#endif

class MnmpSWClient : public blaze::Client {

 public:
  MnmpSWClient(uint32_t* opt_data, uint32_t num_ref_seq, uint64_t* ref_data, uint32_t* mmi_data);

#ifdef DUPLICATE_OUTPUT
  void setup(std::vector<align_input>& inputs, int num_region);
#else
  void setup(std::vector<align_input>& inputs, std::vector<align_output>& outputs, int num_region);
#endif

  void compute();

  int getMaxInputNumItems(int idx);
  int getMaxInputItemLength(int idx);
  int getMaxInputDataWidth(int idx);

 private:
  uint32_t *opt_data_;
  uint32_t  num_ref_seq_;
  uint64_t *ref_data_;
  uint32_t *mmi_data_;

  int           num_region_;
  fpga_inputs  *inputs_data_;
  fpga_outputs *outputs_data_;

};

#endif