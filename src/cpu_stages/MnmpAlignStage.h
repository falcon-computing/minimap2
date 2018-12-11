#ifndef MNMP_FLOW_CPU_ALIGN_H
#define MNMP_FLOW_CPU_ALIGN_H

#include <string>

#include "kflow/Pipeline.h"
#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"
#include "kflow/SourceStage.h"
#include "kflow/SinkStage.h"

#include "MnmpData.h"
#include "MnmpCpuStages.h"


class MinimapAlign : public kestrelFlow::MapStage<ChainsBatch, AlignsBatch, COMPUTE_DEPTH, OUTPUT_DEPTH> {
 public:
  MinimapAlign(int n=1)
  : kestrelFlow::MapStage<ChainsBatch, AlignsBatch, COMPUTE_DEPTH, OUTPUT_DEPTH>(n) {;}

  AlignsBatch compute(ChainsBatch const &i_chainsBatch);
};

#endif
