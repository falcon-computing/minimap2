#ifndef MARKDUP_H
#define MARKDUP_H

#include "kflow/Pipeline.h"
#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"
#include "kflow/SourceStage.h"
#include "kflow/SinkStage.h"

#include "MnmpData.h"
#include "MnmpCpuStages.h"

#include "samblaster.h"
#include <boost/thread/mutex.hpp>

// mapStage version
class MarkDupStage: public kestrelFlow::MapStage<
  BamsBatch, BamsBatch, INPUT_DEPTH, OUTPUT_DEPTH> {
public:
  MarkDupStage(int n = 1, bam_hdr_t * hdr = NULL):kestrelFlow::MapStage<
    BamsBatch, BamsBatch, INPUT_DEPTH, OUTPUT_DEPTH>(n){
      InitializeState(hdr);
      hdr_ = hdr;
    }
  ~MarkDupStage() {
    // deleteState(state);
  } 
  BamsBatch compute(BamsBatch const & input);
private:
  void InitializeState(bam_hdr_t* hdr); 
  state_t* state;
  bam_hdr_t* hdr_;
  boost::mutex mtx_;
};

#endif
