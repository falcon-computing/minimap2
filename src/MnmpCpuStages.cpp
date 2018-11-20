#include "MnmpCpuStages.h"

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

ChainsBatch MinimapChain::compute(SeqsBatch const &i_seqsBatch) {
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started MinimapChain";

  fragExtSOA *l_fragExtSOA;
  if (g_mnmpOpt->flag & MM_F_INDEPEND_SEG) {
    l_fragExtSOA = createFragmentExtensionSOA(i_seqsBatch.m_numSeq);
  }
  else {
    l_fragExtSOA = createFragmentExtensionSOA(i_seqsBatch.m_numFrag);
  }
  for (int l_fr = 0; l_fr < i_seqsBatch.m_numFrag; l_fr++) {
    int l_segOff = i_seqsBatch.m_segOff[l_fr], pe_ori = g_mnmpOpt->pe_ori;
    int qlens[MM_MAX_SEG];
    const char *qseqs[MM_MAX_SEG];
    assert(i_seqsBatch.m_numSeg[l_fr] <= MM_MAX_SEG);
    //if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
    //  DLOG(INFO) << "QR\t" << i_seqsBatch.m_seqs[l_segOffset].name << "\t" << getTid() << "\t" << i_seqsBatch.m_seqs[l_segOffset].l_seq;

    for (int l_sg = 0; l_sg < i_seqsBatch.m_numSeg[l_fr]; ++l_sg) {
      if (i_seqsBatch.m_numSeg[l_fr] == 2 && ((l_sg == 0 && (pe_ori>>1&1)) || (l_sg == 1 && (pe_ori&1))))
        mm_revcomp_bseq(&i_seqsBatch.m_seqs[l_segOff + l_sg]);
      qlens[l_sg] = i_seqsBatch.m_seqs[l_segOff + l_sg].l_seq;
      qseqs[l_sg] = i_seqsBatch.m_seqs[l_segOff + l_sg].seq;
    }

    if (g_mnmpOpt->flag & MM_F_INDEPEND_SEG) {
      for (int l_sg = 0; l_sg < i_seqsBatch.m_numSeg[l_fr]; ++l_sg) {
        int l_repLen, l_fragGap;
        fc_map_frag_chain(g_mnmpOpt, g_minimizer, l_segOff+l_sg, 1, &qlens[l_sg], &qseqs[l_sg], &i_seqsBatch.m_numReg[l_segOff+l_sg], &i_seqsBatch.m_reg[l_segOff+l_sg], i_seqsBatch.m_seqs[l_segOff+l_sg].name, l_fragExtSOA, &l_repLen, &l_fragGap);
        i_seqsBatch.m_repLen[l_segOff + l_sg] = l_repLen;
        i_seqsBatch.m_fragGap[l_segOff + l_sg] = l_fragGap;
      }
    } else {
      int l_repLen, l_fragGap;
      fc_map_frag_chain(g_mnmpOpt, g_minimizer, l_fr, i_seqsBatch.m_numSeg[l_fr], qlens, qseqs, &i_seqsBatch.m_numReg[l_segOff], &i_seqsBatch.m_reg[l_segOff], i_seqsBatch.m_seqs[l_segOff].name, l_fragExtSOA, &l_repLen, &l_fragGap);
      for (int l_sg = 0; l_sg < i_seqsBatch.m_numSeg[l_fr]; ++l_sg) {
        i_seqsBatch.m_repLen[l_segOff + l_sg] = l_repLen;
        i_seqsBatch.m_fragGap[l_segOff + l_sg] = l_fragGap;
      }
    }
  }

  ChainsBatch o_chainsBatch;
  o_chainsBatch.m_batchIdx    = i_seqsBatch.m_batchIdx;
  o_chainsBatch.m_startSeqIdx = i_seqsBatch.m_startSeqIdx;
  o_chainsBatch.m_numSeq      = i_seqsBatch.m_numSeq;
  o_chainsBatch.m_numReg      = i_seqsBatch.m_numReg;
  o_chainsBatch.m_seqs        = i_seqsBatch.m_seqs;
  o_chainsBatch.m_reg         = i_seqsBatch.m_reg;

  o_chainsBatch.m_numFrag     = i_seqsBatch.m_numFrag;
  o_chainsBatch.m_numSeg      = i_seqsBatch.m_numSeg;
  o_chainsBatch.m_segOff      = i_seqsBatch.m_segOff;
  o_chainsBatch.m_repLen      = i_seqsBatch.m_repLen;
  o_chainsBatch.m_fragGap     = i_seqsBatch.m_fragGap;

  o_chainsBatch.m_fragExtSOA  = l_fragExtSOA;

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished MinimapChain";

  return o_chainsBatch;
}

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


void Reorder::compute(int i_workerId) {
  std::map<int, AlignsBatch> l_inorderCache;
  
  std::vector<AlignsBatch> *l_bundle = new std::vector<AlignsBatch>;
  int l_batchCounter = 0;
  int l_bundleCounter = 0;

  while (true) {
    AlignsBatch i_alignsBatch;
    bool l_ready = this->getInput(i_alignsBatch);
    while (!this->isFinal() && !l_ready) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(10));
      l_ready = this->getInput(i_alignsBatch);
    }
    if (!l_ready) {
      break;
    }

    if (!FLAGS_inorder_output) {
      l_bundle->push_back(i_alignsBatch);
      if (l_bundle->size() >= FLAGS_output_size) {
        AlignsBundle o_alignsBundle;
        o_alignsBundle.m_bundleIdx = l_bundleCounter++;
        o_alignsBundle.m_batches = l_bundle;
        this->pushOutput(o_alignsBundle);
        l_bundle = new std::vector<AlignsBatch>;
      }
    }
    else {
      l_inorderCache[i_alignsBatch.m_batchIdx] = i_alignsBatch;
      while (l_inorderCache.find(l_batchCounter) != l_inorderCache.end()) {
        l_bundle->push_back(l_inorderCache[l_batchCounter]);
        l_inorderCache.erase(l_batchCounter);
        l_batchCounter++;
        if (l_bundle->size() >= FLAGS_output_size) {
          AlignsBundle o_alignsBundle;
          o_alignsBundle.m_bundleIdx = l_bundleCounter++;
          o_alignsBundle.m_batches = l_bundle;
          this->pushOutput(o_alignsBundle);
          l_bundle = new std::vector<AlignsBatch>;
        }
      }
    }
  }

  if (l_bundle->size() > 0) {
    AlignsBundle o_alignsBundle;
    o_alignsBundle.m_bundleIdx = l_bundleCounter++;
    o_alignsBundle.m_batches = l_bundle;
    this->pushOutput(o_alignsBundle);
  }
}

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

#if 0
int SeqsWrite::compute(AlignsBatch const &i_alignsBatch) {
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started SeqsWrite";

  const mm_idx_t *mi = g_minimizer;
  void *km = NULL;
  kstring_t l_outputStrBuf = {0, 0, 0};

  // Output File
  std::stringstream l_outputFname;
  l_outputFname << FLAGS_output_dir << "/part-"
                << std::setw(6) << std::setfill('0') << i_alignsBatch.m_batchIdx << ".log";
  std::ofstream l_outputFp;
  samFile      *l_bamOutputFp;
  if (FLAGS_bam) {
    l_bamOutputFp = sam_open(l_outputFname.str().c_str(), "wb0");
  }
  else {
    l_outputFp.open(l_outputFname.str());
  }

  // Write headers
  int l_retCode;
  bam_hdr_t *l_bamHeader = NULL;
  if (FLAGS_bam) {
    std::string l_headerStr = fc_write_sam_hdr(mi, FLAGS_R, VERSION, m_cmdInfo);
    l_bamHeader = sam_hdr_parse(l_headerStr.length(), l_headerStr.c_str());
    l_retCode = sam_hdr_write(l_bamOutputFp, l_bamHeader);
  }
  else
    l_outputFp << fc_write_sam_hdr(mi, FLAGS_R, VERSION, m_cmdInfo);
  
  // Write records
  for (int k = 0; k < i_alignsBatch.m_numFrag; ++k) {
    int seg_st = i_alignsBatch.m_segOff[k];
    int seg_en = i_alignsBatch.m_segOff[k] + i_alignsBatch.m_numSeg[k];
    for (int i = seg_st; i < seg_en; ++i) {
      mm_bseq1_t *t = &i_alignsBatch.m_seqs[i];

      /* NOTE: split prefix is not supported
      if (g_mnmpOpt->split_prefix && p->n_parts == 0) { // then write to temporary files
      	mm_err_fwrite(&i_alignsBatch.m_numReg[i],  sizeof(int), 1, p->fp_split);
      	mm_err_fwrite(&i_alignsBatch.m_repLen[i],  sizeof(int), 1, p->fp_split);
      	mm_err_fwrite(&i_alignsBatch.m_fragGap[i], sizeof(int), 1, p->fp_split);
      	for (j = 0; j < i_alignsBatch.m_numReg[i]; ++j) {
      	  mm_reg1_t *r = &i_alignsBatch.m_reg[i][j];
      	  mm_err_fwrite(r, sizeof(mm_reg1_t), 1, p->fp_split);
      	  if (g_mnmpOpt->flag & MM_F_CIGAR) {
      	    mm_err_fwrite(&r->p->capacity, 4, 1, p->fp_split);
      	    mm_err_fwrite(r->p, r->p->capacity, 4, p->fp_split);
      	  }
      	}
      }
      else */
      if (i_alignsBatch.m_numReg[i] > 0) { // the query has at least one hit
      	for (int j = 0; j < i_alignsBatch.m_numReg[i]; ++j) {
      	  mm_reg1_t *r = &i_alignsBatch.m_reg[i][j];
      	  assert(!r->sam_pri || r->id == r->parent);
      	  if ((g_mnmpOpt->flag & MM_F_NO_PRINT_2ND) && r->id != r->parent)
      	    continue;
          if (FLAGS_bam) {
      	    mm_write_sam2(&l_outputStrBuf, mi, t, i - seg_st, j, i_alignsBatch.m_numSeg[k], &i_alignsBatch.m_numReg[seg_st], (const mm_reg1_t*const*)&i_alignsBatch.m_reg[seg_st], km, g_mnmpOpt->flag);
            bam1_t *l_bam = bam_init1();
            l_retCode = sam_parse1(&l_outputStrBuf, l_bamHeader, l_bam);
            l_retCode = sam_write1(l_bamOutputFp, l_bamHeader, l_bam);
            bam_destroy1(l_bam);
          }
      	  else {
            if (g_mnmpOpt->flag & MM_F_OUT_SAM)
      	      mm_write_sam2(&l_outputStrBuf, mi, t, i - seg_st, j, i_alignsBatch.m_numSeg[k], &i_alignsBatch.m_numReg[seg_st], (const mm_reg1_t*const*)&i_alignsBatch.m_reg[seg_st], km, g_mnmpOpt->flag);
      	    else
      	      mm_write_paf(&l_outputStrBuf, mi, t, r, km, g_mnmpOpt->flag);
      	    //mm_err_puts(l_outputStrBuf.s);
      	    l_outputFp << l_outputStrBuf.s << std::endl;
          }
      	}
      }
      else if (g_mnmpOpt->flag & (MM_F_OUT_SAM|MM_F_PAF_NO_HIT)) { // output an empty hit, if requested
        if (FLAGS_bam) {
      	  mm_write_sam2(&l_outputStrBuf, mi, t, i - seg_st, -1, i_alignsBatch.m_numSeg[k], &i_alignsBatch.m_numReg[seg_st], (const mm_reg1_t*const*)&i_alignsBatch.m_reg[seg_st], km, g_mnmpOpt->flag);
          bam1_t *l_bam = bam_init1();
          l_retCode = sam_parse1(&l_outputStrBuf, l_bamHeader, l_bam);
          l_retCode = sam_write1(l_bamOutputFp, l_bamHeader, l_bam);
          bam_destroy1(l_bam);
        }
      	else {
          if (g_mnmpOpt->flag & MM_F_OUT_SAM)
      	    mm_write_sam2(&l_outputStrBuf, mi, t, i - seg_st, -1, i_alignsBatch.m_numSeg[k], &i_alignsBatch.m_numReg[seg_st], (const mm_reg1_t*const*)&i_alignsBatch.m_reg[seg_st], km, g_mnmpOpt->flag);
      	  else
            mm_write_paf(&l_outputStrBuf, mi, t, 0, 0, g_mnmpOpt->flag);
      	  //mm_err_puts(l_outputStrBuf.s);
      	  l_outputFp << l_outputStrBuf.s << std::endl;
        }
      }
    }
    for (int i = seg_st; i < seg_en; ++i) {
      mm_bseq1_t *l_curSeq = &i_alignsBatch.m_seqs[i];

      for (int j = 0; j < i_alignsBatch.m_numReg[i]; ++j)
        free(i_alignsBatch.m_reg[i][j].p);
      free(i_alignsBatch.m_reg[i]);

      free(l_curSeq->seq);
      free(l_curSeq->name);
      if (l_curSeq->qual)
        free(l_curSeq->qual);
      if (l_curSeq->comment)
        free(l_curSeq->comment);
    }
  }
  free(i_alignsBatch.m_reg);
  free(i_alignsBatch.m_numReg);
  free(i_alignsBatch.m_seqs); // seg_off, n_seg, rep_len and frag_gap were allocated with reg; no memory leak here

  km_destroy(km);
  if (mm_verbose >= 3)
    fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %d sequences\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), i_alignsBatch.m_numSeq);

  if (FLAGS_bam) {
    bam_hdr_destroy(l_bamHeader);
    l_retCode = sam_close(l_bamOutputFp);
  }

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished SeqsWrite";

  return 0;
}
#endif

int SeqsWrite::compute(BamsBatch const &i_bamsBatch) {
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started SeqsWrite";
  static char *l_writeMode[] = { "wb", "wb0", "w" };

  // Output File
  std::stringstream l_outputFname;
  samFile          *l_bamOutputFp;

  l_outputFname << FLAGS_output_dir << "/part-"
                << std::setw(6) << std::setfill('0') << i_bamsBatch.m_batchIdx << ".log";
  if (FLAGS_output_flag < 0 || FLAGS_output_flag > 2)
    FLAGS_output_flag = 1;
  l_bamOutputFp = sam_open(l_outputFname.str().c_str(), l_writeMode[FLAGS_output_flag]);

  // Write headers
  int l_retCode;
  l_retCode = sam_hdr_write(l_bamOutputFp, g_bamHeader);
  
  // Write records
  for (int l_bm = 0; l_bm < i_bamsBatch.m_numBams; l_bm++) {
    l_retCode = sam_write1(l_bamOutputFp, g_bamHeader, i_bamsBatch.m_bams[l_bm]);
    bam_destroy1(i_bamsBatch.m_bams[l_bm]);
  }
  free(i_bamsBatch.m_bams);
  
  l_retCode = sam_close(l_bamOutputFp);

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished SeqsWrite";

  return 0;
}
