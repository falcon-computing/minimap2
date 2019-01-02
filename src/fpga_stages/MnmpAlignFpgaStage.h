#ifndef MNMP_FLOW_FPGA_ALIGN_H
#define MNMP_FLOW_FPGA_ALIGN_H

#include <string>

#include "kflow/Pipeline.h"
#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"
#include "kflow/SourceStage.h"
#include "kflow/SinkStage.h"

#include "MnmpData.h"
#include "MnmpStagesCommon.h"


class MinimapAlignFpga : public kestrelFlow::MapPartitionStage<ChainsBatch, AlignsBatch, COMPUTE_DEPTH, OUTPUT_DEPTH> {
 public:
  MinimapAlignFpga(int n=1)
  : kestrelFlow::MapPartitionStage<ChainsBatch, AlignsBatch, COMPUTE_DEPTH, OUTPUT_DEPTH>(n, false)
  {
    m_num = n;
    m_ttlNumFragsArr = new long[n];
    m_ttlNumFragsOnFPGAArr = new long[n];
  }

  ~MinimapAlignFpga() {
    int l_totalNum = 0, l_fpgaNum = 0;
    for (int i = 0; i < m_num; i++) {
      l_totalNum += m_ttlNumFragsArr[i];
      l_fpgaNum  += m_ttlNumFragsOnFPGAArr[i];
    }
    DLOG(INFO) << "FPGA Alignment Percentage: " << l_fpgaNum << " / " << l_totalNum;
    delete [] m_ttlNumFragsArr;
    delete [] m_ttlNumFragsOnFPGAArr;
  }

  void compute(int i_workerId);
  long *m_ttlNumFragsArr;
  long *m_ttlNumFragsOnFPGAArr;
  int m_num;
};

#endif
