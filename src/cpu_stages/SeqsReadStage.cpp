#include "SeqsReadStage.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>

#include "MnmpGlobal.h"
#include "MnmpOptions.h"
#include "MnmpWrapper.h"
#include "MnmpUtils.h"
#include "minimap.h"
#include "mmpriv.h"
#include "kalloc.h"

#include "htslib/sam.h"


void SeqsRead::compute() {
  // Config reading param
  int l_withQual = (!!(g_mnmpOpt->flag & MM_F_OUT_SAM) && !(g_mnmpOpt->flag & MM_F_NO_QUAL));
  int l_withComment = !!(g_mnmpOpt->flag & MM_F_COPY_COMMENT);
  int l_fragMode = (m_numFp > 1 || !!(g_mnmpOpt->flag & MM_F_FRAG_MODE));
  int l_numProcessed = 0;

  // Open files
  if (m_numFp < 1) {
    return;
  }
  mm_bseq_file_t **l_fp = open_bseqs(m_numFp, (const char **)m_fn);
  if (l_fp == NULL) {
    LOG(ERROR) << "Failed to open input files.";
    return;
  }

  // Read sequences
  int l_batchCounter = 0;
  while (true) {
    DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started SeqsRead";
    SeqsBatch o_seqsBatch;

    if (m_numFp > 1) {
      o_seqsBatch.m_seqs = mm_bseq_read_frag2(m_numFp,
                                              l_fp,
                                              g_mnmpOpt->mini_batch_size,
                                              l_withQual,
                                              l_withComment,
                                             &o_seqsBatch.m_numSeq);
    }
    else {
      o_seqsBatch.m_seqs = mm_bseq_read3(l_fp[0],
                                         g_mnmpOpt->mini_batch_size,
                                         l_withQual,
                                         l_withComment,
                                         l_fragMode,
                                        &o_seqsBatch.m_numSeq);
    }

    if (!o_seqsBatch.m_seqs)
      break;

    o_seqsBatch.m_startSeqIdx = l_numProcessed;
    o_seqsBatch.m_batchIdx    = l_batchCounter++;
    for (int i = 0; i < o_seqsBatch.m_numSeq; i++) {
      o_seqsBatch.m_seqs[i].rid = l_numProcessed++;
    }

    o_seqsBatch.m_numReg = (int*)calloc(5 * o_seqsBatch.m_numSeq, sizeof(int));
    o_seqsBatch.m_segOff = o_seqsBatch.m_numReg + o_seqsBatch.m_numSeq;
    o_seqsBatch.m_numSeg = o_seqsBatch.m_segOff + o_seqsBatch.m_numSeq;
    o_seqsBatch.m_repLen = o_seqsBatch.m_numSeg + o_seqsBatch.m_numSeq;
    o_seqsBatch.m_fragGap = o_seqsBatch.m_repLen + o_seqsBatch.m_numSeq;

    o_seqsBatch.m_reg = (mm_reg1_t**)calloc(o_seqsBatch.m_numSeq, sizeof(mm_reg1_t*));
    //o_seqsBatch.m_buf = (mm_tbuf_t**)calloc(FLAG_t, sizeof(mm_tbuf_t*));

    o_seqsBatch.m_numFrag = 0;
    for (int i = 1, j = 0; i <= o_seqsBatch.m_numSeq; i++) {
      if (i == o_seqsBatch.m_numSeq                                               ||
          !l_fragMode                                                             ||
          !mm_qname_same(o_seqsBatch.m_seqs[i-1].name, o_seqsBatch.m_seqs[i].name)  ) {
        o_seqsBatch.m_numSeg[o_seqsBatch.m_numFrag] = i - j;
        o_seqsBatch.m_segOff[o_seqsBatch.m_numFrag++] = j;
        j = i;
      }
    }

    DLOG_IF(INFO, VLOG_IS_ON(3)) << "Read " << o_seqsBatch.m_numSeq << " sequences";
    DLOG_IF(INFO, VLOG_IS_ON(4)) << "#Frag: " << o_seqsBatch.m_numFrag << ", #Segm: " << o_seqsBatch.m_numSeg[0];
    DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished SeqsRead";

    this->pushOutput(o_seqsBatch);
  }

  // Close files
  for (int l_i = 0; l_i < m_numFp; l_i++)
    mm_bseq_close(l_fp[l_i]);
  free(l_fp);

  return;
}
