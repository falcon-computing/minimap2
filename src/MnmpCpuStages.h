#ifndef MNMP_FLOW_CPU_PIPELINE_H
#define MNMP_FLOW_CPU_PIPELINE_H

#include <string>

#include "kflow/Pipeline.h"
#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"
#include "kflow/SourceStage.h"
#include "kflow/SinkStage.h"

#include "MnmpData.h"


#define INPUT_DEPTH   64
#define OUTPUT_DEPTH  64
#define COMPUTE_DEPTH 64


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

class MinimapOriginMap : public kestrelFlow::MapStage<SeqsBatch, AlignsBatch, INPUT_DEPTH, OUTPUT_DEPTH> {
 public:
  MinimapOriginMap(int n=1)
  : kestrelFlow::MapStage<SeqsBatch, AlignsBatch, INPUT_DEPTH, OUTPUT_DEPTH>(n) {;}

  AlignsBatch compute(SeqsBatch const &i_seqsBatch);
};

class MinimapChain : public kestrelFlow::MapStage<SeqsBatch, ChainsBatch, INPUT_DEPTH, COMPUTE_DEPTH> {
 public:
  MinimapChain(int n=1)
  : kestrelFlow::MapStage<SeqsBatch, ChainsBatch, INPUT_DEPTH, COMPUTE_DEPTH>(n) {;}

  ChainsBatch compute(SeqsBatch const &i_seqsBatch);
};

class MinimapAlign : public kestrelFlow::MapStage<ChainsBatch, AlignsBatch, COMPUTE_DEPTH, OUTPUT_DEPTH> {
 public:
  MinimapAlign(int n=1)
  : kestrelFlow::MapStage<ChainsBatch, AlignsBatch, COMPUTE_DEPTH, OUTPUT_DEPTH>(n) {;}

  AlignsBatch compute(ChainsBatch const &i_chainsBatch);
};

class Reorder : public kestrelFlow::MapPartitionStage<AlignsBatch, AlignsBundle, COMPUTE_DEPTH, COMPUTE_DEPTH> {
 public:
  Reorder(int n=1)
  : kestrelFlow::MapPartitionStage<AlignsBatch, AlignsBundle, COMPUTE_DEPTH, COMPUTE_DEPTH>(n, false) {;}

  void compute(int i_workerId); 
};

class CoordSort : public kestrelFlow::MapStage<AlignsBundle, BamsBatch, COMPUTE_DEPTH, COMPUTE_DEPTH> {
 public:
  CoordSort(int n=1)
  : kestrelFlow::MapStage<AlignsBundle, BamsBatch, COMPUTE_DEPTH, COMPUTE_DEPTH>(n) {;}

  BamsBatch compute(AlignsBundle const &i_alignsBundle);
};

#if 0
class SeqsWrite : public kestrelFlow::MapStage<AlignsBatch, int, OUTPUT_DEPTH, 0> {
 public:
  SeqsWrite(int n=1, std::string i_cmdInfo = "")
  : kestrelFlow::MapStage<AlignsBatch, int, OUTPUT_DEPTH, 0>(n),
    m_cmdInfo(i_cmdInfo)
  {;}

  int compute(AlignsBatch const &i_alignsBatch);
 private:
  std::string m_cmdInfo;
};
#endif

class SeqsWrite : public kestrelFlow::MapStage<BamsBatch, int, OUTPUT_DEPTH, 0> {
 public:
  SeqsWrite(int n=1)
  : kestrelFlow::MapStage<BamsBatch, int, OUTPUT_DEPTH, 0>(n) {;}

  int compute(BamsBatch const &i_bamsBatch);
};

#endif
