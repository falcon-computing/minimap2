#ifndef MNMP_FLOW_CPU_KTOB_H
#define MNMP_FLOW_CPU_KTOB_H

#include <string>

#include "kflow/Pipeline.h"
#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"
#include "kflow/SourceStage.h"
#include "kflow/SinkStage.h"

#include "MnmpData.h"
#include "MnmpCpuStages.h"

class KseqsToBseqs : public kestrelFlow::MapStage<KseqsBatch, SeqsBatch, INPUT_DEPTH, COMPUTE_DEPTH> {
  public:
    KseqsToBseqs(int i_numFp, int n=1)
    : kestrelFlow::MapStage<KseqsBatch, SeqsBatch, INPUT_DEPTH, COMPUTE_DEPTH>(n),
    m_numFp(i_numFp) {;}

    SeqsBatch compute(KseqsBatch const &i_kseqsBatch);
  private:
    int  m_numFp;
};

#endif