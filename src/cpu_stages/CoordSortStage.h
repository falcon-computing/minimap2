#ifndef MNMP_FLOW_CPU_COORDSORT_H
#define MNMP_FLOW_CPU_COORDSORT_H

#include <string>

#include "kflow/Pipeline.h"
#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"
#include "kflow/SourceStage.h"
#include "kflow/SinkStage.h"

#include "MnmpData.h"
#include "MnmpCpuStages.h"


class CoordSort : public kestrelFlow::MapStage<AlignsBundle, BamsBatch, COMPUTE_DEPTH, COMPUTE_DEPTH> {
 public:
  CoordSort(int n=1)
  : kestrelFlow::MapStage<AlignsBundle, BamsBatch, COMPUTE_DEPTH, COMPUTE_DEPTH>(n) {;}

  BamsBatch compute(AlignsBundle const &i_alignsBundle);
};

#endif
