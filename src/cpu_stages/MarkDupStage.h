#ifndef MARKDUP
#define MARKDUP
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
  AlignsBundle, BamsBatch, INPUT_DEPTH, OUTPUT_DEPTH> {
public:
  MarkDupStage(int n = 1, bam_hdr_t* head = NULL):kestrelFlow::MapStage<
    AlignsBundle, BamsBatch, INPUT_DEPTH, OUTPUT_DEPTH>(n), head_(head){
      InitializeState(head_);
    }
  ~MarkDupStage() {
    // deleteState(state);
  } 
  BamsBatch compute(AlignsBundle const & input);
private:
  void InitializeState(bam_hdr_t* head); 
  state_t* state_;
  bam_hdr_t* head_;
  boost::mutex mtx_;
};
#endif
