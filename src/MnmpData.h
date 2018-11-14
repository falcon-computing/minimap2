#ifndef MNMP_FLOW_DATA_H
#define MNMP_FLOW_DATA_H

#include <string>
#include "bseq.h"
#include "minimap.h"


struct fragExtSOA {
  mm_reg1_t **m_regs0Arr;
  int        *m_numRegs0;

  uint32_t  *m_hash;
  mm128_t  **m_anchorArr;
};

fragExtSOA *createFragmentExtensionSOA(int);
void deleteFragmentExtensionSOA(fragExtSOA *);

struct SeqsBatch {
  int         m_batchIdx;
  int         m_startSeqIdx;
  int         m_numSeq; // n_seq
  mm_bseq1_t *m_seqs;

  int  m_numFrag;
  int *m_numReg;
  int *m_segOff;
  int *m_numSeg;
  int *m_repLen;
  int *m_fragGap;

  mm_reg1_t **m_reg;

  //std::string name_tag = "SeqsBatch"; 
};

struct ChainsBatch {
  int         m_batchIdx;
  int         m_startSeqIdx;
  int         m_numSeq; // n_seq
  mm_bseq1_t *m_seqs;

  int  m_numFrag;
  int *m_numReg;
  int *m_segOff;
  int *m_numSeg;
  int *m_repLen;
  int *m_fragGap;

  mm_reg1_t **m_reg;

  fragExtSOA *m_fragExtSOA;

  //std::string name_tag = "ChainsBatch";
};


struct AlignsBatch {
  int         m_batchIdx;
  int         m_startSeqIdx;
  int         m_numSeq; // n_seq
  mm_bseq1_t *m_seqs;

  int  m_numFrag;
  int *m_numReg;
  int *m_segOff;
  int *m_numSeg;
  int *m_repLen;
  int *m_fragGap;

  mm_reg1_t **m_reg;

  //std::string name_tag = "AlignsBatch";
};




#endif
