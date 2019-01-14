#ifndef MNMP_FLOW_CPU_READ_H
#define MNMP_FLOW_CPU_READ_H

#include <string>

#include "kflow/Pipeline.h"
#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"
#include "kflow/SourceStage.h"
#include "kflow/SinkStage.h"

#include "MnmpData.h"
#include "MnmpCpuStages.h"
#include "MnmpUtils.h"

class SeqsRead : public kestrelFlow::SourceStage<SeqsBatch, INPUT_DEPTH> {
 public:
  SeqsRead(int i_numFp, char **i_fn) : kestrelFlow::SourceStage<SeqsBatch, INPUT_DEPTH>(),
                                       m_numFp(i_numFp),
                                       m_fn(i_fn)
  {;}
  void compute();
 private:
  int    m_numFp;
  char **m_fn;
};

#endif
