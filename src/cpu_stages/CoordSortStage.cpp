#include "CoordSortStage.h"

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

static bool bam1_lt(const bam1_t *a, const bam1_t *b) {
  return ((uint64_t)a->core.tid<<32|(a->core.pos+1)<<1|bam_is_rev(a))
       < ((uint64_t)b->core.tid<<32|(b->core.pos+1)<<1|bam_is_rev(b));
}

BamsBatch CoordSort::compute(AlignsBundle const &i_alignsBundle) {
  bam_hdr_t *l_bamHeader = g_bamHeader;
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started CoordSort";

  int l_numBamsEst = 0;
  for (int l_ba = 0; l_ba < i_alignsBundle.m_batches->size(); l_ba++) {
    AlignsBatch &l_alignsBatch = (*i_alignsBundle.m_batches)[l_ba];
    for (int l_sq = 0; l_sq < l_alignsBatch.m_numSeq; l_sq++) {
      l_numBamsEst += (l_alignsBatch.m_numReg[l_sq] > 0) ? l_alignsBatch.m_numReg[l_sq]
                      : ((g_mnmpOpt->flag & (MM_F_OUT_SAM|MM_F_PAF_NO_HIT)) ? 1 : 0);
    }
  }
  l_numBamsEst = l_numBamsEst + (int)(0.2*l_numBamsEst);

  // Convert into bams
  kstring_t l_samStrBuf = {0, 0, 0};
  bam1_t **l_bamsArr = (bam1_t**)malloc(l_numBamsEst*sizeof(bam1_t*));
  int l_numBams = 0;
  for (int l_ba = 0; l_ba < i_alignsBundle.m_batches->size(); l_ba++) {
    AlignsBatch &l_alignsBatch = (*i_alignsBundle.m_batches)[l_ba];
    for (int l_fr = 0; l_fr < l_alignsBatch.m_numFrag; l_fr++) {
      int l_segOffset = l_alignsBatch.m_segOff[l_fr];
      int l_numSegs   = l_alignsBatch.m_numSeg[l_fr];

      int *l_numRegsArr = &l_alignsBatch.m_numReg[l_segOffset];
      mm_reg1_t const *const * l_regsArr = &l_alignsBatch.m_reg[l_segOffset];

      // Convert a fragment[pair] to bams
      int l_retCode;
      for (int l_sg = l_segOffset; l_sg < l_segOffset+l_numSegs; l_sg++) {
        mm_bseq1_t *l_seq = &l_alignsBatch.m_seqs[l_sg];

        if (l_alignsBatch.m_numReg[l_sg] > 0) { // the query has at least one hit
          for (int l_rg = 0; l_rg < l_alignsBatch.m_numReg[l_sg]; l_rg++) {
            mm_reg1_t *l_reg = &l_alignsBatch.m_reg[l_sg][l_rg];
            assert(!l_reg->sam_pri || l_reg->id == l_reg->parent);
            if ((g_mnmpOpt->flag & MM_F_NO_PRINT_2ND) && l_reg->id != l_reg->parent)
              continue;

            mm_write_sam2(&l_samStrBuf, g_minimizer, l_seq, l_sg-l_segOffset, l_rg, l_numSegs, l_numRegsArr, l_regsArr, NULL, g_mnmpOpt->flag);
            if (l_numBams >= l_numBamsEst) {
              l_numBamsEst += (int)(0.2*l_numBamsEst);
              l_bamsArr = (bam1_t**)realloc(l_bamsArr, l_numBamsEst*sizeof(bam1_t*));
            }
            bam1_t *l_bam = bam_init1();
            l_retCode = sam_parse1(&l_samStrBuf, l_bamHeader, l_bam);
            l_bamsArr[l_numBams++] = l_bam;
          }
        }
        else if (g_mnmpOpt->flag & (MM_F_OUT_SAM|MM_F_PAF_NO_HIT)) { // output an empty hit, if requested
          mm_write_sam2(&l_samStrBuf, g_minimizer, l_seq, l_sg-l_segOffset, -1, l_numSegs, l_numRegsArr, l_regsArr, NULL, g_mnmpOpt->flag);
          if (l_numBams >= l_numBamsEst) {
            l_numBamsEst += (int)(0.2*l_numBamsEst);
            l_bamsArr = (bam1_t**)realloc(l_bamsArr, l_numBamsEst*sizeof(bam1_t*));
          }
          bam1_t *l_bam = bam_init1();
          l_retCode = sam_parse1(&l_samStrBuf, l_bamHeader, l_bam);
          l_bamsArr[l_numBams++] = l_bam;
        }
      }

      // Free a fragment[pair]; TODO: retain for PAF output
      for (int l_sg = l_segOffset; l_sg < l_segOffset+l_numSegs; l_sg++) {
        mm_bseq1_t *l_curSeq = &l_alignsBatch.m_seqs[l_sg];

        for (int l_rg = 0; l_rg < l_alignsBatch.m_numReg[l_sg]; l_rg++)
          free(l_alignsBatch.m_reg[l_sg][l_rg].p);
        free(l_alignsBatch.m_reg[l_sg]);

        free(l_curSeq->seq);
        free(l_curSeq->name);
        if (l_curSeq->qual)
          free(l_curSeq->qual);
        if (l_curSeq->comment)
          free(l_curSeq->comment);
      }
    } // close loop over fragments

    // Free a batch of alignment
    free(l_alignsBatch.m_reg);
    free(l_alignsBatch.m_numReg);
    free(l_alignsBatch.m_seqs); // seg_off, n_seg, rep_len and frag_gap were allocated with reg; no memory leak here
  } // close loop over batches
  delete i_alignsBundle.m_batches;
  l_bamsArr = (bam1_t**)realloc(l_bamsArr, l_numBams*sizeof(bam1_t*));
  free(l_samStrBuf.s);


  // Sort
  if (FLAGS_sort) {
    std::sort(l_bamsArr, l_bamsArr+l_numBams, bam1_lt);
  }

  BamsBatch o_bamsBatch;
  o_bamsBatch.m_batchIdx = i_alignsBundle.m_bundleIdx;
  o_bamsBatch.m_bams     = l_bamsArr;
  o_bamsBatch.m_numBams  = l_numBams;
  
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished CoordSort";

  return o_bamsBatch;
}
