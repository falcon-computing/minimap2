#include "MnmpChainStage.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <thread>

#include <sched.h>
#include <numa.h>

#include "MnmpGlobal.h"
#include "MnmpOptions.h"
#include "MnmpWrapper.h"
#include "MnmpUtils.h"
#include "minimap.h"
#include "mmpriv.h"
#include "kalloc.h"

#include "htslib/sam.h"



ChainsBatch MinimapChain::compute(SeqsBatch const &i_seqsBatch) {
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started MinimapChain";

  // A simple HARDCODED routine to determine NUMA node
  // TODO: parse output of `lscpu -p` to get NUMA node info
  thread_local const unsigned int C_NUM_TOTAL_CPU  = std::thread::hardware_concurrency();
  thread_local const int          C_NUM_TOTAL_NODE = (numa_available()==-1) ? 1 : numa_num_configured_nodes();
  
  int l_cpuId  = sched_getcpu();
  int l_coreId = l_cpuId % (C_NUM_TOTAL_CPU/2); // assume USE_SMT:2
  int l_nodeId = l_coreId / (C_NUM_TOTAL_CPU/2/C_NUM_TOTAL_NODE);

  mm_idx_t *l_minimizer;
  if (FLAGS_use_numa)
    l_minimizer = g_minimizerNumaList[l_nodeId];
  else
    l_minimizer = g_minimizer;

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
        fc_map_frag_chain(g_mnmpOpt, l_minimizer, l_segOff+l_sg, 1, &qlens[l_sg], &qseqs[l_sg], &i_seqsBatch.m_numReg[l_segOff+l_sg], &i_seqsBatch.m_reg[l_segOff+l_sg], i_seqsBatch.m_seqs[l_segOff+l_sg].name, l_fragExtSOA, &l_repLen, &l_fragGap);
        i_seqsBatch.m_repLen[l_segOff + l_sg] = l_repLen;
        i_seqsBatch.m_fragGap[l_segOff + l_sg] = l_fragGap;
      }
    } else {
      int l_repLen, l_fragGap;
      fc_map_frag_chain(g_mnmpOpt, l_minimizer, l_fr, i_seqsBatch.m_numSeg[l_fr], qlens, qseqs, &i_seqsBatch.m_numReg[l_segOff], &i_seqsBatch.m_reg[l_segOff], i_seqsBatch.m_seqs[l_segOff].name, l_fragExtSOA, &l_repLen, &l_fragGap);
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
