#ifndef MNMP_FLOW_CPU_REORDER_H
#define MNMP_FLOW_CPU_REORDER_H

#include <string>

#include "kflow/Pipeline.h"
#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"
#include "kflow/SourceStage.h"
#include "kflow/SinkStage.h"

#include "MnmpData.h"
#include "MnmpCpuStages.h"


class Reorder : public kestrelFlow::MapPartitionStage<AlignsBatch, AlignsBundle, COMPUTE_DEPTH, COMPUTE_DEPTH> {
 public:
  Reorder(int n=1)
  : kestrelFlow::MapPartitionStage<AlignsBatch, AlignsBundle, COMPUTE_DEPTH, COMPUTE_DEPTH>(n, false) {;}

  void compute(int i_workerId); 
};

#endif
