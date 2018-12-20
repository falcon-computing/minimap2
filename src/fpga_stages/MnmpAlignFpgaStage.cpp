#include "MnmpAlignFpgaStage.h"

#include <vector>
#include <queue>
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

#include "./MnmpAlignFPGAUtils/host_interface/host_types.h"
// #ifndef LOCAL_BLAZE
#include "./MnmpAlignFPGAUtils/host_interface/mmsw_orig.h"
#include "./MnmpAlignFPGAUtils/host_interface/mmsw_baseline.h"
#include "./MnmpAlignFPGAUtils/host_interface/mmsw_fpga.h"
// #else
#include "./MnmpAlignFPGAUtils/blaze_interface/MnmpSWWorker.h"
#include "./MnmpAlignFPGAUtils/blaze_interface/MnmpSWClient.h"
// #endif


#include "htslib/sam.h"

bool packInput(int const           i_fragIdx,
               int const           i_numSegs,
               mm_bseq1_t const *  i_segSeqs,
               mm_seg_t         *  i_segChains,
               int              *  io_numRegsOfSgmArr,
               mm_reg1_t        ** io_regsOfSgmArr,
               std::vector<align_input> &o_fpgaChunkInput,
               std::vector<align_input> &o_fpgaChunkInputNext )
{
  extern unsigned char seq_nt4_table[256];

  thread_local std::queue<align_input> tl_packBuffer;
  thread_local int                     tl_queryLens[MM_MAX_SEG];
  thread_local char const *            tl_querySeqs[MM_MAX_SEG];
  thread_local uint8_t                 tl_encodedQseq[2][MAX_SEQ_LENGTH];

  for (int l_sg = 0; l_sg < i_numSegs; l_sg++) {
    if (i_segSeqs[l_sg].l_seq > MAX_SEQ_LENGTH)
      return false;
    
    tl_queryLens[l_sg] = i_segSeqs[l_sg].l_seq;
    tl_querySeqs[l_sg] = i_segSeqs[l_sg].seq;
  }
  
  for (int l_sg = 0; i_numSegs > 1 && l_sg < i_numSegs; l_sg++) {
    /* the number of  regions of a segment in the fragment, aka a sequence in the pair */
    int         l_numRegs = io_numRegsOfSgmArr[l_sg];
    /* pointer to regions of a segment in the fragment, aka a sequence in the pair */
    mm_reg1_t * l_regs    = io_regsOfSgmArr[l_sg];

    // update mm_reg1_t::parent (original scope: mm_map_frag/fc_map_frag_align)
    mm_set_parent(NULL,
                  g_mnmpOpt->mask_level,
                  l_numRegs,
                  l_regs,
                  g_mnmpOpt->a * 2 + g_mnmpOpt->b,
                  g_mnmpOpt->flag&MM_F_HARD_MLEVEL);
    
    // alias anchors (original scope: align_regs)
    // mm128_t *l_anchors = (i_numSegs > 1) ? i_segChains[l_sg].a : i_segChains->a;
    mm128_t *l_anchors = i_segChains[l_sg].a;
    // get squeezed number of anchors (original scope: mm_align_skeleton, after query seqs encoding)
    int32_t l_numAnchors = mm_squeeze_a(NULL, l_numRegs, l_regs, l_anchors);

    if (l_numAnchors > MAX_ANKOR_NUM) {
      std::queue<align_input> l_emptyQ;
      std::swap(tl_packBuffer, l_emptyQ);
      return false;
    }

    // encode the query sequence (original scope: mm_align_skeleton, before regions squeezing)
    int  const         l_qlen = tl_queryLens[l_sg];
    char const * const l_qseq = tl_querySeqs[l_sg];
    for (int l_bs = 0; l_bs < l_qlen; l_bs++) {
      tl_encodedQseq[0][l_bs] = seq_nt4_table[(uint8_t)l_qseq[l_bs]];
      tl_encodedQseq[1][tl_queryLens[l_sg]-l_bs-1] = (tl_encodedQseq[0][l_bs] < 4) ? 3-tl_encodedQseq[0][l_bs] : 4;
    }
    
    for (int l_rg = 0; l_rg < l_numRegs; l_rg++) {
      align_input l_newAlignInput;
      l_newAlignInput.qlen = l_qlen;
      for (int l_bs = 0; l_bs < l_qlen; l_bs++) {
        l_newAlignInput.qseq[0][l_bs] = tl_encodedQseq[0][l_bs];
        l_newAlignInput.qseq[1][l_bs] = tl_encodedQseq[1][l_bs];
      }

      l_newAlignInput.n_a = l_numAnchors;
      memcpy(l_newAlignInput.a, l_anchors, l_numAnchors*sizeof(mm128_t));

      /* leave region[1] as blank */
      l_newAlignInput.region[0].orig = l_regs[l_rg];
      /* leave dp_score[1] as blank */
      // l_newAlignInput.dp_score[0] = l_regs[l_rg].p->dp_score;

      l_newAlignInput.data_size = sizeof(align_input) - 16 * MAX_ANKOR_NUM + 16 * l_numAnchors;

      tl_packBuffer.push(l_newAlignInput);
    }
  }

  if (tl_packBuffer.size() > MAX_BATCH_SIZE) {
    std::queue<align_input> l_emptyQ;
    std::swap(tl_packBuffer, l_emptyQ);
    return false;
  }
  
  if (tl_packBuffer.size() + o_fpgaChunkInput.size() <= MAX_BATCH_SIZE) {
    while (!tl_packBuffer.empty()) {
      o_fpgaChunkInput.push_back(tl_packBuffer.front());
      tl_packBuffer.pop();
    }
  }
  else {
    while (!tl_packBuffer.empty()) {
      o_fpgaChunkInputNext.push_back(tl_packBuffer.front());
      tl_packBuffer.pop();
    }
  }
  
  return true;
}

bool unpackOutput(int const           i_fragIdx,
                  int const           i_numSegs,
                  mm_bseq1_t       *  i_segSeqs,
                  mm_seg_t         *  i_segChains,
                  int const           i_repLen,
                  int const           i_fragGap,
                  int              *  io_numRegsOfSgmArr,
                  mm_reg1_t        ** io_regsOfSgmArr,
                  align_input      *  i_alignInputs,
                  align_output     *  i_alignOutputs)
{
  thread_local int const is_sr = !!(g_mnmpOpt->flag & MM_F_SR);
  thread_local int const pe_ori = g_mnmpOpt->pe_ori;

  thread_local int tl_queryLensPair[2];

  thread_local uint8_t *tl_encQseqPtr[MM_MAX_SEG];

  ksw_extz_t l_ez;

  int l_resOff = 0;
  for (int l_sg = 0; l_sg < i_numSegs; l_sg++) {
    memset(&l_ez, 0, sizeof(ksw_extz_t));

    bool l_needAlign = false;
    for (int l_rg = 0; l_rg < io_numRegsOfSgmArr[l_sg]; l_rg++) {
      align_input  * l_alignInput;
      align_output * l_alignOutput;
      mm_reg1_t * l_reg = &( io_regsOfSgmArr[l_sg][l_rg] );
      mm_reg1_t   l_reg2;

      // Copy output from FPGA results
      if (!l_needAlign) {
        l_alignInput  = &i_alignInputs[l_resOff];
        l_alignOutput = &i_alignOutputs[l_resOff++];
        // l_reg = &( io_regsOfSgmArr[l_sg][l_rg] );
        *l_reg = l_alignOutput->region[0].orig;

        size_t l_cigarSize = l_alignOutput->p.n_cigar*sizeof(uint32_t);
        l_reg->p = (mm_extra_t*)realloc(l_reg->p, sizeof(mm_extra_t)+l_cigarSize);
        *(l_reg->p) = l_alignOutput->p;
        memcpy(l_reg->p->cigar, l_alignOutput->cigar, l_cigarSize); 

        l_reg2 = l_alignOutput->region[1].orig;
        l_reg2.p = NULL;
      }
      // Do alignment on CPU
      else {
        assert( l_rg != 0 );
        assert( l_resOff > 0 );
        tl_encQseqPtr[0] = &(i_alignInputs[l_resOff-1].qseq[0][0]);
        tl_encQseqPtr[1] = &(i_alignInputs[l_resOff-1].qseq[1][0]);
        mm_align1(NULL,
                  g_mnmpOpt,
                  g_minimizer,
                  i_alignInputs[l_resOff-1].qlen,
                  tl_encQseqPtr,
                  l_reg,
                 &l_reg2,
                  i_alignInputs[l_resOff-1].n_a,
                  i_alignInputs[l_resOff-1].a,
                 &l_ez,
                  g_mnmpOpt->flag);
      }

      if (g_mnmpOpt->flag & MM_F_SPLICE)
        l_reg->p->trans_strand = g_mnmpOpt->flag&MM_F_SPLICE_FOR ? 1 : 2;
      
      // Insert new region, marked to be processed next
      if (l_reg2.cnt > 0) {
        io_regsOfSgmArr[l_sg] = mm_insert_reg(&l_reg2,
                                               l_rg,
                                              &io_numRegsOfSgmArr[l_sg],
                                               io_regsOfSgmArr[l_sg]);
        l_needAlign = true; // do alignment in the next iteration
      }
      else {
        l_needAlign = false;
      }
      
      // Realias region pointer as reallocation in *mm_insert_reg*
      l_reg = &( io_regsOfSgmArr[l_sg][l_rg] );

      // Insert new region for inversion, marked to be skipped
      if (l_rg > 0 && l_reg->split_inv) {
        tl_encQseqPtr[0] = &(l_alignInput->qseq[0][0]);
        tl_encQseqPtr[1] = &(l_alignInput->qseq[1][0]);
  			if (mm_align1_inv(NULL,
                          g_mnmpOpt,
                          g_minimizer,
                          l_alignInput->qlen,
                          tl_encQseqPtr,
                          l_reg-1,
                          l_reg,
                         &l_reg2,
                         &l_ez)               ) {
  				io_regsOfSgmArr[l_sg] = mm_insert_reg(&l_reg2,
                                                 l_rg,
                                                &io_numRegsOfSgmArr[l_sg],
                                                 io_regsOfSgmArr[l_sg]);
  				++l_rg; // skip the inserted INV alignment
  			}
  		}
    } // close loop over regions

    if (l_ez.cigar != NULL)
      kfree(NULL, l_ez.cigar);
    assert( l_resOff > 0 || io_numRegsOfSgmArr[l_sg] == 0 );
    mm_filter_regs(g_mnmpOpt,
                   (io_numRegsOfSgmArr[l_sg] == 0) ? 0 : i_alignInputs[l_resOff-1].qlen,
                  &io_numRegsOfSgmArr[l_sg],
                   io_regsOfSgmArr[l_sg]);
    mm_hit_sort(NULL, &io_numRegsOfSgmArr[l_sg], io_regsOfSgmArr[l_sg]);

    // Don't choose primary mapping(s)
    if (!(g_mnmpOpt->flag & MM_F_ALL_CHAINS)) {
  		mm_set_parent(NULL,
                    g_mnmpOpt->mask_level,
                    io_numRegsOfSgmArr[l_sg],
                    io_regsOfSgmArr[l_sg],
                    g_mnmpOpt->a * 2 + g_mnmpOpt->b,
                    g_mnmpOpt->flag&MM_F_HARD_MLEVEL);
  		mm_select_sub(NULL,
                    g_mnmpOpt->pri_ratio,
                    g_minimizer->k*2,
                    g_mnmpOpt->best_n,
                   &io_numRegsOfSgmArr[l_sg],
                    io_regsOfSgmArr[l_sg]);
  		mm_set_sam_pri(io_numRegsOfSgmArr[l_sg], io_regsOfSgmArr[l_sg]);
  	}

    mm_set_mapq(NULL, io_numRegsOfSgmArr[l_sg], io_regsOfSgmArr[l_sg], g_mnmpOpt->min_chain_score, g_mnmpOpt->a, i_repLen, is_sr);
  } // close loop over segments

  mm_seg_free(NULL, i_numSegs, i_segChains);
  if (i_numSegs == 2               &&
      g_mnmpOpt->pe_ori >= 0       &&
      (g_mnmpOpt->flag&MM_F_CIGAR)   ) {
    tl_queryLensPair[0] = i_segSeqs[0].l_seq;
    tl_queryLensPair[1] = i_segSeqs[1].l_seq;
    mm_pair(NULL, i_fragGap, g_mnmpOpt->pe_bonus, g_mnmpOpt->a * 2 + g_mnmpOpt->b, g_mnmpOpt->a, tl_queryLensPair, io_numRegsOfSgmArr, io_regsOfSgmArr);
  }

  // Flip the query strand and coordinate to the original read strand
  for (int l_sg = 0; l_sg < i_numSegs; ++l_sg) {
    if (i_numSegs == 2                    &&
        ((l_sg == 0 && (pe_ori>>1&1)) ||
         (l_sg == 1 && (pe_ori&1))      )   ) {
      mm_revcomp_bseq(&i_segSeqs[l_sg]);
      for (int l_rg = 0; l_rg < io_numRegsOfSgmArr[l_sg]; l_rg++) {
        mm_reg1_t *r = &io_regsOfSgmArr[l_sg][l_rg];
        int t = r->qs;
        r->qs = i_segSeqs[l_sg].l_seq - r->qe;
        r->qe = i_segSeqs[l_sg].l_seq - t;
        r->rev = !r->rev;
      }
    }
  }
}

static void alignmentOnCpu(int l_fr, ChainsBatch const &i_chainsBatch, fragExtSOA *l_fragExtSOA) {
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

#define ALIGN_WORKER_MSG(_worker_idx_) ("[AlignFPGA "+std::to_string(_worker_idx_)+ "] ")
void MinimapAlignFpga::compute(int i_workerId) {
  long tl_totalNumFrags = 0;
  long tl_totalNumFragsOnFPGA = 0;
  DLOG(INFO) << ALIGN_WORKER_MSG(i_workerId) << "Initialized FPGA worker for Alignment";

  std::string l_btsm = "/curr/jyqiu/workspace/minimap2_release/mmsw_kernel_xilinx_vcu1525_dynamic_5_1_1211.xclbin";
#ifndef LOCAL_BLAZE
  mmsw_fpga_init(g_mnmpOpt, &l_btsm[0], g_minimizer);
  DLOG(INFO) << "Initialized FPGA";
#else
  mm_mapopt_t_union* l_fpgaOpt = (mm_mapopt_t_union*)new mm_mapopt_t_union;
  mm_idx_fpga*       l_fpgaIdx = (mm_idx_fpga*)new mm_idx_fpga; 

  l_fpgaOpt->orig = *g_mnmpOpt;
  ref_conversion(g_minimizer, l_fpgaIdx);

  MnmpSWClient l_client(l_fpgaOpt->serial, l_fpgaIdx->n_seq, l_fpgaIdx->seqs, l_fpgaIdx->S);
#endif

  const int l_fpgaChunkSize = 2048;

  while (true) {
    ChainsBatch i_chainsBatch;
    bool l_ready = this->getInput(i_chainsBatch);
    while (!this->isFinal() && !l_ready) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(5));
      l_ready = this->getInput(i_chainsBatch);
    }
    if (!l_ready) {
      break;
    }

    DLOG_IF(INFO, VLOG_IS_ON(1)) << ALIGN_WORKER_MSG(i_workerId) << "Started MinimapAlign";

    // for frIdx = 0 -> NUM_FRAGMENT
    //     NUM_SEGMENT := NumSegsArr[frIdx]
    //     sgmOff := SgmOffArr[frIdx]
    //     repLen := RepLenArr[sgmOff]
    //     fragGap := FragGapArr[sgmOff]
    //
    //     for sgmIdx = 0 -> NUM_SEGMENT
    //         seqIdx := sgmOff + sgmIdx
    //

    // MULTI-WAY BUFFER PATTERN:
    // KERNEL      |-K1--|      |--K1------|      |--K1------|      |--K1--|
    // IN-OUT    INA INB OUTA INC OUTB INA OUTC INB OUTA INC OUTB INA OUTC OUTA
    // KERNEL          |--K2------|      |--K2------|      |--K2-------|

    int        * const l_numSegsArr    = i_chainsBatch.m_numSeg;
    int        * const l_segOffArr     = i_chainsBatch.m_segOff;
    int        * const l_repLenArr     = i_chainsBatch.m_repLen;
    int        * const l_fragGapArr    = i_chainsBatch.m_fragGap;
    mm_seg_t  ** const l_segChainsArr  = i_chainsBatch.m_fragExtSOA->m_segChainsArr;
    int        * const l_numRegsArr    = i_chainsBatch.m_numReg;
    mm_reg1_t ** const l_regsArr       = i_chainsBatch.m_reg;

    std::vector<int>          l_fragIds[3];
    std::vector<align_input>  l_fpgaChunkInput[3];
    std::vector<align_output> l_fpgaChunkOutput[3];

    int l_startFragIdx = 0;
    int l_endFragIdx = 0;
    int l_prevFragIdx = -1;
    int l_bufferIdx = 0;
    int l_startPacking = 1;
    for (int l_fr = 0; l_fr < i_chainsBatch.m_numFrag; l_fr++) {
      if (l_startPacking) {
        DLOG_IF(INFO, VLOG_IS_ON(3)) << ALIGN_WORKER_MSG(i_workerId) << "Pack input into buffer " << l_bufferIdx;
        l_startPacking = 0;
      }
      int         const l_numSegs  = l_numSegsArr[l_fr];
      int         const l_segOff   = l_segOffArr[l_fr];
      int         const l_repLen   = l_repLenArr[l_fr];
      int         const l_fragGap  = l_fragGapArr[l_fr];
      mm_seg_t  * const l_segChains= l_segChainsArr[l_fr];
      int       *       l_numRegsOfSgmArr = &l_numRegsArr[l_segOff];
      mm_reg1_t **      l_regsOfSgmArr    = &l_regsArr[l_segOff];

      // Pack input
      bool l_isPacked = packInput(l_fr,
                                  l_numSegs,
                                 &i_chainsBatch.m_seqs[l_segOff],
                                  l_segChains,
                                  l_numRegsOfSgmArr,
                                  l_regsOfSgmArr,
                                  l_fpgaChunkInput[l_bufferIdx],
                                  l_fpgaChunkInput[(l_bufferIdx+1)%3]);
      if (l_isPacked) {
        if (l_fpgaChunkInput[(l_bufferIdx+1)%3].size() > 0)
          l_fragIds[(l_bufferIdx+1)%3].push_back(l_fr);
        else
          l_fragIds[l_bufferIdx].push_back(l_fr);
      }
      
      // Continue packing
      if (l_fpgaChunkInput[l_bufferIdx].size() < MAX_BATCH_SIZE &&
          l_fpgaChunkInput[(l_bufferIdx+1)%3].size() == 0         ) {
        continue;
      }

      // Invoke kernel
      DLOG_IF(INFO, VLOG_IS_ON(3)) << ALIGN_WORKER_MSG(i_workerId) << "Call kernel for buffer " << l_bufferIdx;
      // mmsw_orig_compute(g_mnmpOpt, g_minimizer, l_fpgaChunkInput[l_bufferIdx], l_fpgaChunkOutput[l_bufferIdx]);
      // double c, d;
      // mmsw_baseline_compute(g_mnmpOpt, g_minimizer, l_fpgaChunkInput[l_bufferIdx], l_fpgaChunkOutput[l_bufferIdx], c, d);
#ifndef LOCAL_BLAZE
      mmsw_fpga_compute(&l_btsm[0], g_mnmpOpt, g_minimizer, l_fpgaChunkInput[l_bufferIdx], l_fpgaChunkOutput[l_bufferIdx]);
#else
      MnmpSWWorker l_worker(&l_client, l_fpgaChunkInput[l_bufferIdx], l_fpgaChunkInput[l_bufferIdx].size());
      l_worker.run();
      l_worker.getOutput(l_fpgaChunkOutput[l_bufferIdx]);
#endif

      // Switch buffer
      l_bufferIdx = (l_bufferIdx+1)%3;

      // Unpack output
      int l_resOff = 0;
      if (l_fragIds[(l_bufferIdx+2)%3].size() > 0)
        DLOG_IF(INFO, VLOG_IS_ON(3)) << ALIGN_WORKER_MSG(i_workerId) << "Unpack output to buffer " << (l_bufferIdx+2)%3;
      for (int l_fr2 : l_fragIds[(l_bufferIdx+2)%3]) {
        for (int l_fr3 = l_prevFragIdx+1; l_fr3 < l_fr2; l_fr3++) {
          alignmentOnCpu(l_fr3, i_chainsBatch, i_chainsBatch.m_fragExtSOA);
        }
        l_prevFragIdx = l_fr2;

        int         const l_numSegs2  = l_numSegsArr[l_fr2];
        int         const l_segOff2   = l_segOffArr[l_fr2];
        int         const l_repLen2   = l_repLenArr[l_fr2];
        int         const l_fragGap2  = l_fragGapArr[l_fr2];
        mm_seg_t  * const l_segChains2= l_segChainsArr[l_fr2];
        int       *       l_numRegsOfSgmArr2 = &l_numRegsArr[l_segOff2];
        mm_reg1_t **      l_regsOfSgmArr2    = &l_regsArr[l_segOff2];
        align_input  * l_alignInputs = &l_fpgaChunkInput[(l_bufferIdx+2)%3][l_resOff];
        align_output * l_alignOutputs = &l_fpgaChunkOutput[(l_bufferIdx+2)%3][l_resOff];

        for (int l_sg = 0; l_sg < l_numSegs2; l_sg++)
          l_resOff += l_numRegsOfSgmArr2[l_sg];

        unpackOutput(l_fr2,
                     l_numSegs2,
                    &i_chainsBatch.m_seqs[l_segOff2],
                     l_segChains2,
                     l_repLen2,
                     l_fragGap2,
                     l_numRegsOfSgmArr2,
                     l_regsOfSgmArr2,
                     l_alignInputs,
                     l_alignOutputs);
        tl_totalNumFragsOnFPGA++;
      }
      l_fragIds[(l_bufferIdx+2)%3].clear();
      l_fpgaChunkInput[(l_bufferIdx+2)%3].clear();
      l_fpgaChunkOutput[(l_bufferIdx+2)%3].clear();

      // // Switch buffer
      // l_bufferIdx = (l_bufferIdx+1)%3;
      l_startPacking = 1;
    }

    // Peeling: Invoke kernel if last input buffer is not empty
    DLOG_IF(INFO, VLOG_IS_ON(3)) << ALIGN_WORKER_MSG(i_workerId) << "Call kernel for buffer " << l_bufferIdx;
    if (l_fpgaChunkInput[l_bufferIdx].size() > 0) {
      // Do alignment on FPGA for a inputs number large enough 
      if (l_fpgaChunkInput[l_bufferIdx].size() >= 256) {
        // mmsw_orig_compute(g_mnmpOpt, g_minimizer, l_fpgaChunkInput[l_bufferIdx], l_fpgaChunkOutput[l_bufferIdx]);
        // double a, b;
        // mmsw_baseline_compute(g_mnmpOpt, g_minimizer, l_fpgaChunkInput[l_bufferIdx], l_fpgaChunkOutput[l_bufferIdx], a, b);
#ifndef LOCAL_BLAZE
      mmsw_fpga_compute(&l_btsm[0], g_mnmpOpt, g_minimizer, l_fpgaChunkInput[l_bufferIdx], l_fpgaChunkOutput[l_bufferIdx]);
#else
      MnmpSWWorker l_worker(&l_client, l_fpgaChunkInput[l_bufferIdx], l_fpgaChunkInput[l_bufferIdx].size());
      l_worker.run();
      l_worker.getOutput(l_fpgaChunkOutput[l_bufferIdx]);
#endif
      }
      else {
        l_fragIds[l_bufferIdx].clear();
        l_fpgaChunkInput[l_bufferIdx].clear();
      }
    }

    // Peeling: Unpack output for last two kernel launches
    l_bufferIdx = (l_bufferIdx+2)%3;
    for (int l_k = 0; l_k < 2; l_k++, l_bufferIdx = (l_bufferIdx+1)%3) {
      if (l_fragIds[l_bufferIdx].size() > 0)
        DLOG_IF(INFO, VLOG_IS_ON(3)) << ALIGN_WORKER_MSG(i_workerId) << "Unpack output to buffer " << l_bufferIdx;
      int l_resOff = 0;
      for (int l_fr2 : l_fragIds[l_bufferIdx]) {
        for (int l_fr3 = l_prevFragIdx+1; l_fr3 < l_fr2; l_fr3++) {
          alignmentOnCpu(l_fr3, i_chainsBatch, i_chainsBatch.m_fragExtSOA);
        }
        l_prevFragIdx = l_fr2;

        int         const l_numSegs2  = l_numSegsArr[l_fr2];
        int         const l_segOff2   = l_segOffArr[l_fr2];
        int         const l_repLen2   = l_repLenArr[l_fr2];
        int         const l_fragGap2  = l_fragGapArr[l_fr2];
        mm_seg_t  * const l_segChains2= l_segChainsArr[l_fr2];
        int       *       l_numRegsOfSgmArr2 = &l_numRegsArr[l_segOff2];
        mm_reg1_t **      l_regsOfSgmArr2    = &l_regsArr[l_segOff2];
        align_input  * l_alignInputs = &l_fpgaChunkInput[l_bufferIdx][l_resOff];
        align_output * l_alignOutputs = &l_fpgaChunkOutput[l_bufferIdx][l_resOff];

        for (int l_sg = 0; l_sg < l_numSegs2; l_sg++)
          l_resOff += l_numRegsOfSgmArr2[l_sg];

        unpackOutput(l_fr2,
                     l_numSegs2,
                    &i_chainsBatch.m_seqs[l_segOff2],
                     l_segChains2,
                     l_repLen2,
                     l_fragGap2,
                     l_numRegsOfSgmArr2,
                     l_regsOfSgmArr2,
                     l_alignInputs,
                     l_alignOutputs);
        tl_totalNumFragsOnFPGA++;
      }
      l_fragIds[l_bufferIdx].clear();
      l_fpgaChunkInput[l_bufferIdx].clear();
      l_fpgaChunkOutput[l_bufferIdx].clear();
    }

    // Peeling: Do alignment on CPU for tailing unpacked fragments
    for (int l_fr3 = l_prevFragIdx+1; l_fr3 < i_chainsBatch.m_numFrag; l_fr3++)
      alignmentOnCpu(l_fr3, i_chainsBatch, i_chainsBatch.m_fragExtSOA);

    deleteFragmentExtensionSOA(i_chainsBatch.m_fragExtSOA);

    tl_totalNumFrags += i_chainsBatch.m_numFrag;

#if 0
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

        int  l_numSegs =  i_chainsBatch.m_numSeg[l_fr];
        int *l_numRegs = &i_chainsBatch.m_numReg[l_segOff];
        mm_reg1_t **regs = &i_chainsBatch.m_reg[l_segOff];
        {
          uint32_t l_hash = l_fragExtSOA->m_hash[l_fr];
          mm_reg1_t *l_regs0Arr = l_fragExtSOA->m_regs0Arr[l_fr];
          int l_numRegs0 = l_fragExtSOA->m_numRegs0[l_fr];
          mm128_t *l_anchorArr = l_fragExtSOA->m_anchorArr[l_fr];

          int is_sr = !!(g_mnmpOpt->flag & MM_F_SR);

          if (l_numSegs == 1) { // uni-segment
            l_regs0Arr = align_regs(g_mnmpOpt, g_minimizer, NULL, qlens[0], qseqs[0], &l_numRegs0, l_regs0Arr, l_anchorArr);
            mm_set_mapq(NULL, l_numRegs0, l_regs0Arr, g_mnmpOpt->g_minimizern_chain_score, g_mnmpOpt->a, l_repLen, is_sr);
            l_numRegs[0] = l_numRegs0, regs[0] = l_regs0Arr;
          } else { // multi-segment
            mm_seg_t *seg;
            seg = mm_seg_gen(NULL, l_hash, l_numSegs, qlens, l_numRegs0, l_regs0Arr, l_numRegs, regs, l_anchorArr); // split fragment chain to separate segment chains        
            free(l_regs0Arr);
            for (int i = 0; i < l_numSegs; ++i) {
              mm_set_parent(NULL, g_mnmpOpt->mask_level, l_numRegs[i], regs[i], g_mnmpOpt->a * 2 + g_mnmpOpt->b, g_mnmpOpt->flag&MM_F_HARD_MLEVEL); // update mm_reg1_t::parent 
              regs[i] = align_regs(g_mnmpOpt, g_minimizer, NULL, qlens[i], qseqs[i], &l_numRegs[i], regs[i], seg[i].a);
              mm_set_mapq(NULL, l_numRegs[i], regs[i], g_mnmpOpt->g_minimizern_chain_score, g_mnmpOpt->a, l_repLen, is_sr);
            }           
            mm_seg_free(NULL, l_numSegs, seg);
            if (l_numSegs == 2 && g_mnmpOpt->pe_ori >= 0 && (g_mnmpOpt->flag&MM_F_CIGAR))
              mm_pair(NULL, l_fragGap, g_mnmpOpt->pe_bonus, g_mnmpOpt->a * 2 + g_mnmpOpt->b, g_mnmpOpt->a, qlens, l_numRegs, regs); // pairing
          }                   

          kfree(NULL, l_anchorArr);
        }
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
#endif

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

    DLOG_IF(INFO, VLOG_IS_ON(1)) << ALIGN_WORKER_MSG(i_workerId) << "Finished MinimapAlign";

    this->pushOutput(o_alignsBatch); 

  }

#ifndef LOCAL_BLAZE
  mmsw_fpga_destroy();
#else
  if (!l_fpgaOpt) 
    delete l_fpgaOpt;
  if (!l_fpgaIdx) 
    delete l_fpgaIdx;
#endif

  DLOG(INFO) << ALIGN_WORKER_MSG(i_workerId) << "Retired FPGA worker for Alignment";

  m_ttlNumFragsArr[i_workerId] = tl_totalNumFrags;
  m_ttlNumFragsOnFPGAArr[i_workerId] = tl_totalNumFragsOnFPGA;
  return;
}
