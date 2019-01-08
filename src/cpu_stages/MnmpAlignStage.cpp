#include "MnmpAlignStage.h"

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


AlignsBatch MinimapAlign::compute(ChainsBatch const &i_chainsBatch) {
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started MinimapAlign";

  fragExtSOA *l_fragExtSOA = i_chainsBatch.m_fragExtSOA;
  for (int l_fr = 0; l_fr < i_chainsBatch.m_numFrag; l_fr++) {
    int l_segOff = i_chainsBatch.m_segOff[l_fr], pe_ori = g_mnmpOpt->pe_ori;
    int qlens[MM_MAX_SEG];
    const char *qseqs[MM_MAX_SEG];
    for (int l_sg = 0; l_sg < i_chainsBatch.m_numSeg[l_fr]; ++l_sg) {
      qlens[l_sg] = i_chainsBatch.m_seqs[l_segOff + l_sg].l_seq;
      qseqs[l_sg] = i_chainsBatch.m_seqs[l_segOff + l_sg].seq;
    }

    if (g_mnmpOpt->flag & MM_F_INDEPEND_SEG) {
      for (int l_sg = 0; l_sg < i_chainsBatch.m_numSeg[l_fr]; ++l_sg) {
        int l_repLen = i_chainsBatch.m_repLen[l_segOff + l_sg];
        int l_fragGap = i_chainsBatch.m_fragGap[l_segOff + l_sg];
        fc_map_frag_align(g_mnmpOpt, g_minimizer, l_segOff+l_sg, 1, &qlens[l_sg], &qseqs[l_sg], &i_chainsBatch.m_numReg[l_segOff+l_sg], &i_chainsBatch.m_reg[l_segOff+l_sg], i_chainsBatch.m_seqs[l_segOff+l_sg].name, l_repLen, l_fragGap, l_fragExtSOA);
      }
    } else {
      int l_repLen = i_chainsBatch.m_repLen[l_segOff];
      int l_fragGap = i_chainsBatch.m_fragGap[l_segOff];
      fc_map_frag_align(g_mnmpOpt, g_minimizer, l_fr, i_chainsBatch.m_numSeg[l_fr], qlens, qseqs, &i_chainsBatch.m_numReg[l_segOff], &i_chainsBatch.m_reg[l_segOff], i_chainsBatch.m_seqs[l_segOff].name, l_repLen, l_fragGap, l_fragExtSOA);
    }

    for (int l_sg = 0; l_sg < i_chainsBatch.m_numSeg[l_fr]; ++l_sg) {// flip the query strand and coordinate to the original read strand
      if (i_chainsBatch.m_numSeg[l_fr] == 2 && ((l_sg == 0 && (pe_ori>>1&1)) || (l_sg == 1 && (pe_ori&1)))) {
        int k, t;
        mm_revcomp_bseq(&i_chainsBatch.m_seqs[l_segOff + l_sg]);
        for (k = 0; k < i_chainsBatch.m_numReg[l_segOff + l_sg]; ++k) {
          mm_reg1_t *r = &i_chainsBatch.m_reg[l_segOff + l_sg][k];
          t = r->qs;
          r->qs = qlens[l_sg] - r->qe;
          r->qe = qlens[l_sg] - t;
          r->rev = !r->rev;
        }
      }
    }
  }
  deleteFragmentExtensionSOA(l_fragExtSOA);

  AlignsBatch o_alignsBatch;
  o_alignsBatch.m_batchIdx    = i_chainsBatch.m_batchIdx;
  o_alignsBatch.m_startSeqIdx = i_chainsBatch.m_startSeqIdx;
  o_alignsBatch.m_numSeq      = i_chainsBatch.m_numSeq;
  o_alignsBatch.m_numReg      = i_chainsBatch.m_numReg;
  o_alignsBatch.m_seqs        = i_chainsBatch.m_seqs;
  o_alignsBatch.m_reg         = i_chainsBatch.m_reg;

  o_alignsBatch.m_numFrag     = i_chainsBatch.m_numFrag;
  o_alignsBatch.m_numSeg      = i_chainsBatch.m_numSeg;
  o_alignsBatch.m_segOff      = i_chainsBatch.m_segOff;
  o_alignsBatch.m_repLen      = i_chainsBatch.m_repLen;
  o_alignsBatch.m_fragGap     = i_chainsBatch.m_fragGap;

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished MinimapAlign";

  return o_alignsBatch; 
}
