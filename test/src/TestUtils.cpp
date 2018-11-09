#include "TestUtils.h"

#include "MnmpGlobal.h"
#include "MnmpWrapper.h"

void getSeqsBatch(SeqsBatch &io_seqsBatch) {
  int l_withQual = (!!(g_mnmpOpt->flag & MM_F_OUT_SAM) && !(g_mnmpOpt->flag & MM_F_NO_QUAL));
  int l_withComment = !!(g_mnmpOpt->flag & MM_F_COPY_COMMENT);
  int l_fragMode = (g_numFp > 1 || !!(g_mnmpOpt->flag & MM_F_FRAG_MODE));

  mm_bseq_file_t **l_fp = open_bseqs(g_numFp, (const char **)g_fn);
  if (g_numFp > 1) {
    io_seqsBatch.m_seqs = mm_bseq_read_frag2(g_numFp,
                                            l_fp,
                                            g_mnmpOpt->mini_batch_size,
                                            l_withQual,
                                            l_withComment,
                                           &io_seqsBatch.m_numSeq);
  }
  else {
    io_seqsBatch.m_seqs = mm_bseq_read3(l_fp[0],
                                       g_mnmpOpt->mini_batch_size,
                                       l_withQual,
                                       l_withComment,
                                       l_fragMode,
                                      &io_seqsBatch.m_numSeq);
  }

  for (int i = 1024; i < io_seqsBatch.m_numSeq; i++) {
    mm_bseq1_t *l_curSeq = &io_seqsBatch.m_seqs[i];
    free(l_curSeq->seq);
    free(l_curSeq->name);
    if (l_curSeq->qual)
      free(l_curSeq->qual);
    if (l_curSeq->comment)
      free(l_curSeq->comment); 
  }
  io_seqsBatch.m_numSeq = std::min(1024, io_seqsBatch.m_numSeq);

  io_seqsBatch.m_startSeqIdx = 0;
  io_seqsBatch.m_batchIdx    = 0;
  for (int i = 0; i < io_seqsBatch.m_numSeq; i++) {
    io_seqsBatch.m_seqs[i].rid = i;
  }

  io_seqsBatch.m_numReg = (int*)calloc(5 * io_seqsBatch.m_numSeq, sizeof(int));
  io_seqsBatch.m_segOff = io_seqsBatch.m_numReg + io_seqsBatch.m_numSeq;
  io_seqsBatch.m_numSeg = io_seqsBatch.m_segOff + io_seqsBatch.m_numSeq;
  io_seqsBatch.m_repLen = io_seqsBatch.m_numSeg + io_seqsBatch.m_numSeq;
  io_seqsBatch.m_fragGap = io_seqsBatch.m_repLen + io_seqsBatch.m_numSeq;

  io_seqsBatch.m_reg = (mm_reg1_t**)calloc(io_seqsBatch.m_numSeq, sizeof(mm_reg1_t*));

  io_seqsBatch.m_numFrag = 0;
  for (int i = 1, j = 0; i <= io_seqsBatch.m_numSeq; i++) {
    if (i == io_seqsBatch.m_numSeq                                                ||
        !l_fragMode                                                               ||
        !mm_qname_same(io_seqsBatch.m_seqs[i-1].name, io_seqsBatch.m_seqs[i].name)  ) {
      io_seqsBatch.m_numSeg[io_seqsBatch.m_numFrag] = i - j;
      io_seqsBatch.m_segOff[io_seqsBatch.m_numFrag++] = j;
      j = i;
    }
  }
}

void releaseAlignsBatch(AlignsBatch &io_alignsBatch) {
  for (int i = 0; i < io_alignsBatch.m_numSeq; i++) {
    mm_bseq1_t *l_curSeq = &io_alignsBatch.m_seqs[i];

    for (int j = 0; j < io_alignsBatch.m_numReg[i]; ++j)
      free(io_alignsBatch.m_reg[i][j].p);
    free(io_alignsBatch.m_reg[i]);

    free(l_curSeq->seq);
    free(l_curSeq->name);
    if (l_curSeq->qual)
      free(l_curSeq->qual);
    if (l_curSeq->comment)
      free(l_curSeq->comment); 
  }

  free(io_alignsBatch.m_seqs);
  free(io_alignsBatch.m_numReg);
  free(io_alignsBatch.m_reg);
}
