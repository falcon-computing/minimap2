#ifndef MNMP_FLOW_CPU_WRITE_H
#define MNMP_FLOW_CPU_WRITE_H

#include <string>

#include "kflow/Pipeline.h"
#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"
#include "kflow/SourceStage.h"
#include "kflow/SinkStage.h"

#include "MnmpData.h"
#include "MnmpCpuStages.h"


class SeqsWrite : public kestrelFlow::MapStage<BamsBatch, int, OUTPUT_DEPTH, 0> {
 public:
  SeqsWrite(int n=1)
  : kestrelFlow::MapStage<BamsBatch, int, OUTPUT_DEPTH, 0>(n) {;}

  int compute(BamsBatch const &i_bamsBatch);
};

#if 0
class SeqsWrite : public kestrelFlow::MapStage<AlignsBatch, int, OUTPUT_DEPTH, 0> {
 public:
  SeqsWrite(int n=1, std::string i_cmdInfo = "")
  : kestrelFlow::MapStage<AlignsBatch, int, OUTPUT_DEPTH, 0>(n),
    m_cmdInfo(i_cmdInfo)
  {;}

  int compute(AlignsBatch const &i_alignsBatch);
 private:
  std::string m_cmdInfo;
};
#endif

#endif
