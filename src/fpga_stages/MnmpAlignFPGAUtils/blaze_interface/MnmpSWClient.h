#ifndef MMSWCLIENT_H
#define MMSWCLIENT_H
#include <stdexcept>

#include "blaze/Client.h"
#include "ksight/tools.h"
#include "../MnmpAlignFPGAUtils/host_interface/host_types.h"
#include "../MnmpAlignFPGAUtils/host_interface/mmsw_fpga.h"

class MnmpSWClient : public blaze::Client {

 public:
  MnmpSWClient(uint32_t* opt_data, uint32_t num_ref_seq, uint64_t* ref_data, uint32_t* mmi_data);

  void setup(std::vector<align_input>& inputs, int num_region);

  void compute();

  int getMaxInputNumItems(int idx);
  int getMaxInputItemLength(int idx);
  int getMaxInputDataWidth(int idx);

 private:
  uint32_t num_ref_seq_;

  int num_region_;

};

#endif