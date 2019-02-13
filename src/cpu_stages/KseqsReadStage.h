#ifndef MNMP_FLOW_CPU_KREAD_H
#define MNMP_FLOW_CPU_KREAD_H

#include <string>

#include "kflow/Pipeline.h"
#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"
#include "kflow/SourceStage.h"
#include "kflow/SinkStage.h"

#include "MnmpData.h"
#include "MnmpCpuStages.h"

class KseqsRead : public kestrelFlow::SourceStage<KseqsBatch, INPUT_DEPTH> {
 public:
  KseqsRead(int i_numFp, char **i_fn);
  void compute();
  ~KseqsRead();
 private:
  int  m_numFp;
  char **m_fn;
};

#endif