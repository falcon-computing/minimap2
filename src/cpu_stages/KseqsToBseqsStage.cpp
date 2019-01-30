#include "KseqsToBseqsStage.h"

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

// copied from bseq.c
static inline char *kstrdup(const kstring_t *s) {
	char *t;
	t = (char*)malloc(s->l + 1);
	memcpy(t, s->s, s->l + 1);
	return t;
}

// copied from bseq.c
static inline void kseq2bseq(kseq_t *ks, mm_bseq1_t *s, int with_qual, int with_comment) {
	int i;
	s->name = kstrdup(&ks->name);
	s->seq = kstrdup(&ks->seq);
	for (i = 0; i < (int)ks->seq.l; ++i) // convert U to T
		if (s->seq[i] == 'u' || s->seq[i] == 'U')
			--s->seq[i];
	s->qual = with_qual && ks->qual.l? kstrdup(&ks->qual) : 0;
	s->comment = with_comment && ks->comment.l? kstrdup(&ks->comment) : 0;
	s->l_seq = ks->seq.l;
}

SeqsBatch KseqsToBseqs::compute(KseqsBatch const &i_kseqsBatch) {
  DLOG(INFO) << "Started KseqsToBseqsStage";
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "number of seqs in kseqsBatch " << i_kseqsBatch.m_numSeq;
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "m_startSeqIdx of kseqsBatch " << i_kseqsBatch.m_startSeqIdx;
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "m_batchIdx of kseqsBatch " << i_kseqsBatch.m_batchIdx;

  // Config reading param
  int l_withQual = (!!(g_mnmpOpt->flag & MM_F_OUT_SAM) && !(g_mnmpOpt->flag & MM_F_NO_QUAL));
  int l_withComment = !!(g_mnmpOpt->flag & MM_F_COPY_COMMENT);
  int l_fragMode = (m_numFp > 1 || !!(g_mnmpOpt->flag & MM_F_FRAG_MODE));

  SeqsBatch o_seqsBatch;
  mm_bseq1_t * m_seqs = (mm_bseq1_t *)malloc(i_kseqsBatch.m_numSeq * sizeof(mm_bseq1_t));
  o_seqsBatch.m_seqs = m_seqs;

  for (int i = 0; i < i_kseqsBatch.m_numSeq; i ++) {
    kseq2bseq(i_kseqsBatch.kseqs[i], m_seqs + i, l_withQual, l_withComment);
    free(i_kseqsBatch.kseqs[i]->name.s);
    free(i_kseqsBatch.kseqs[i]->comment.s);
    free(i_kseqsBatch.kseqs[i]->seq.s);
    free(i_kseqsBatch.kseqs[i]->qual.s);
  }
  if (i_kseqsBatch.m_numSeq > 0) {
    ks_destroy(i_kseqsBatch.kseqs[0]->f);
  }
  for (int i = 0; i < i_kseqsBatch.m_numSeq; i++) {
    free(i_kseqsBatch.kseqs[i]);
  }
  free(i_kseqsBatch.kseqs);

  o_seqsBatch.m_startSeqIdx = i_kseqsBatch.m_startSeqIdx;
  o_seqsBatch.m_batchIdx = i_kseqsBatch.m_batchIdx;
  o_seqsBatch.m_numSeq = i_kseqsBatch.m_numSeq;
  for (int i = 0; i < o_seqsBatch.m_numSeq; i ++) {
    o_seqsBatch.m_seqs[i].rid = i + o_seqsBatch.m_startSeqIdx;
  }

  o_seqsBatch.m_numReg = (int*)calloc(5 * o_seqsBatch.m_numSeq, sizeof(int));
  o_seqsBatch.m_segOff = o_seqsBatch.m_numReg + o_seqsBatch.m_numSeq;
  o_seqsBatch.m_numSeg = o_seqsBatch.m_segOff + o_seqsBatch.m_numSeq;
  o_seqsBatch.m_repLen = o_seqsBatch.m_numSeg + o_seqsBatch.m_numSeq;
  o_seqsBatch.m_fragGap = o_seqsBatch.m_repLen + o_seqsBatch.m_numSeq;

  o_seqsBatch.m_reg = (mm_reg1_t**)calloc(o_seqsBatch.m_numSeq, sizeof(mm_reg1_t*));

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

  DLOG(INFO) << "Finished KseqsToBseqsStage";

  return o_seqsBatch;
}