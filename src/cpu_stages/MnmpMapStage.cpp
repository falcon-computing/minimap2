#include "MnmpMapStage.h"

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


AlignsBatch MinimapOriginMap::compute(SeqsBatch const &i_seqsBatch) {
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started MnmpOriMap";

  mm_tbuf_t *b = mm_tbuf_init();
  for (int l_fr = 0; l_fr < i_seqsBatch.m_numFrag; l_fr++) {
    int qlens[MM_MAX_SEG], off = i_seqsBatch.m_segOff[l_fr], pe_ori = g_mnmpOpt->pe_ori;
    const char *qseqs[MM_MAX_SEG];
    assert(i_seqsBatch.m_numSeg[l_fr] <= MM_MAX_SEG);
    if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
      DLOG(INFO) << "QR\t" << i_seqsBatch.m_seqs[off].name << "\t" << getTid() << "\t" << i_seqsBatch.m_seqs[off].l_seq;
    for (int l_sg = 0; l_sg < i_seqsBatch.m_numSeg[l_fr]; ++l_sg) {
      if (i_seqsBatch.m_numSeg[l_fr] == 2 && ((l_sg == 0 && (pe_ori>>1&1)) || (l_sg == 1 && (pe_ori&1))))
        mm_revcomp_bseq(&i_seqsBatch.m_seqs[off + l_sg]);
      qlens[l_sg] = i_seqsBatch.m_seqs[off + l_sg].l_seq;
      qseqs[l_sg] = i_seqsBatch.m_seqs[off + l_sg].seq;
    }
    if (g_mnmpOpt->flag & MM_F_INDEPEND_SEG) {
      for (int l_sg = 0; l_sg < i_seqsBatch.m_numSeg[l_fr]; ++l_sg) {
        mm_map_frag(g_minimizer, 1, &qlens[l_sg], &qseqs[l_sg], &i_seqsBatch.m_numReg[off+l_sg], &i_seqsBatch.m_reg[off+l_sg], b, g_mnmpOpt, i_seqsBatch.m_seqs[off+l_sg].name);
        i_seqsBatch.m_repLen[off + l_sg] = b->rep_len;
        i_seqsBatch.m_fragGap[off + l_sg] = b->frag_gap;
      }
    } else {
      mm_map_frag(g_minimizer, i_seqsBatch.m_numSeg[l_fr], qlens, qseqs, &i_seqsBatch.m_numReg[off], &i_seqsBatch.m_reg[off], b, g_mnmpOpt, i_seqsBatch.m_seqs[off].name);
      for (int l_sg = 0; l_sg < i_seqsBatch.m_numSeg[l_fr]; ++l_sg) {
        i_seqsBatch.m_repLen[off + l_sg] = b->rep_len;
        i_seqsBatch.m_fragGap[off + l_sg] = b->frag_gap;
      }
    }
    for (int l_sg = 0; l_sg < i_seqsBatch.m_numSeg[l_fr]; ++l_sg) {// flip the query strand and coordinate to the original read strand
      if (i_seqsBatch.m_numSeg[l_fr] == 2 && ((l_sg == 0 && (pe_ori>>1&1)) || (l_sg == 1 && (pe_ori&1)))) {
        int k, t;
        mm_revcomp_bseq(&i_seqsBatch.m_seqs[off + l_sg]);
        for (k = 0; k < i_seqsBatch.m_numReg[off + l_sg]; ++k) {
          mm_reg1_t *r = &i_seqsBatch.m_reg[off + l_sg][k];
          t = r->qs;
          r->qs = qlens[l_sg] - r->qe;
          r->qe = qlens[l_sg] - t;
          r->rev = !r->rev;
        }
      }
    }
  }
  mm_tbuf_destroy(b);


  AlignsBatch o_alignsBatch;
  o_alignsBatch.m_batchIdx    = i_seqsBatch.m_batchIdx;
  o_alignsBatch.m_startSeqIdx = i_seqsBatch.m_startSeqIdx;
  o_alignsBatch.m_numSeq      = i_seqsBatch.m_numSeq;
  o_alignsBatch.m_numReg      = i_seqsBatch.m_numReg;
  o_alignsBatch.m_seqs        = i_seqsBatch.m_seqs;
  o_alignsBatch.m_reg         = i_seqsBatch.m_reg;

  o_alignsBatch.m_numFrag     = i_seqsBatch.m_numFrag;
  o_alignsBatch.m_numSeg      = i_seqsBatch.m_numSeg;
  o_alignsBatch.m_segOff      = i_seqsBatch.m_segOff;
  o_alignsBatch.m_repLen      = i_seqsBatch.m_repLen;
  o_alignsBatch.m_fragGap     = i_seqsBatch.m_fragGap;

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished MnmpOriMap";

  return o_alignsBatch; 
}
