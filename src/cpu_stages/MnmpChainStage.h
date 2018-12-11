#ifndef MNMP_FLOW_CPU_CHAIN_H
#define MNMP_FLOW_CPU_CHAIN_H

#include <string>

#include "kflow/Pipeline.h"
#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"
#include "kflow/SourceStage.h"
#include "kflow/SinkStage.h"

#include "MnmpData.h"
#include "MnmpCpuStages.h"


class MinimapChain : public kestrelFlow::MapStage<SeqsBatch, ChainsBatch, INPUT_DEPTH, COMPUTE_DEPTH> {
 public:
  MinimapChain(int n=1)
  : kestrelFlow::MapStage<SeqsBatch, ChainsBatch, INPUT_DEPTH, COMPUTE_DEPTH>(n) {;}

  ChainsBatch compute(SeqsBatch const &i_seqsBatch);
};

#endif
