#ifndef MNMP_FLOW_CPU_MAP_H
#define MNMP_FLOW_CPU_MAP_H

#include <string>

#include "kflow/Pipeline.h"
#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"
#include "kflow/SourceStage.h"
#include "kflow/SinkStage.h"

#include "MnmpData.h"
#include "MnmpCpuStages.h"


class MinimapOriginMap : public kestrelFlow::MapStage<SeqsBatch, AlignsBatch, INPUT_DEPTH, OUTPUT_DEPTH> {
 public:
  MinimapOriginMap(int n=1)
  : kestrelFlow::MapStage<SeqsBatch, AlignsBatch, INPUT_DEPTH, OUTPUT_DEPTH>(n) {;}

  AlignsBatch compute(SeqsBatch const &i_seqsBatch);
};

#endif
