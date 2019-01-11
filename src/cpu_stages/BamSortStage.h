#ifndef BAMSORTSTAGE_H
#define BAMSORTSTAGE_H

#include "MnmpData.h"
#include "MnmpCpuStages.h"

class BamSortStage
: public kestrelFlow::MapStage<
      BamRecord, BamRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>
{
 public:
  BamSortStage(int n=1): kestrelFlow::MapStage<
      BamRecord, BamRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>(n) {;} 

  BamRecord compute(BamRecord const & input);
};
#endif
