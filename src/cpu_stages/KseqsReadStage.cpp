#include "KseqsReadStage.h"

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

#define CHECK_PAIR_THRES 1000000

// copied from bseq.c
struct mm_bseq_file_s {
	gzFile fp;
	kseq_t *ks;
	mm_bseq1_t s;
};

void KseqsRead::compute() {

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

  // a buffer kseq for mm_bseq_read3
  kseq_t * tmp_kseq = NULL;

  while (true) {
    DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started KseqsRead";
    KseqsBatch o_kseqsBatch;
    int kseqs_size = 256; // shouldn't be set less than 1
    kseq_t ** kseqs = (kseq_t **)calloc(kseqs_size, sizeof(kseq_t *));
    if (kseqs == NULL) {
      throw("can not allocate kseqs array.\n");
    }
    int kseqs_i = 0;
    int64_t size = 0;

    if (m_numFp > 1) {
      // read as in function mm_bseq_read_frag2
      while(1) {
        int n_read = 0;
        for (int i = 0; i < m_numFp; i++) {
          kseq_t * tmp = kseq_init(l_fp[i]->fp);
          if (kseq_read(tmp) >= 0) {
            n_read++;
            if (kseqs_i >= kseqs_size) {
              kseqs_size *= 2;
              kseqs = (kseq_t **)realloc(kseqs, kseqs_size * sizeof(kseq_t *));
            }
            kseqs[kseqs_i] = tmp;
            kseqs_i ++;
            size += tmp->seq.l;
          }
        }
        if (n_read < m_numFp) {
          if (n_read > 0) {
              DLOG_IF(INFO, VLOG_IS_ON(2)) << "query files have different number of records; extra records skipped.";
          }
          break;
        }
        DLOG_IF(INFO, VLOG_IS_ON(2)) << "size " << size << " g_mnmpOpt->mini_batch_size " << g_mnmpOpt->mini_batch_size;
        if (size >= g_mnmpOpt->mini_batch_size) {
          break;
        }
      }
    }
    else {
      // read as in function mm_bseq_read3
      if (tmp_kseq != NULL) {
        kseqs[0] = tmp_kseq; // should always has 1 space in the beginning
        size += tmp_kseq->seq.l;
        tmp_kseq = NULL;
      }
      while (1) {
        kseq_t * tmp = kseq_init(l_fp[0]->fp);
        if (kseq_read(tmp) < 0) {
          break;
        }
        if (kseqs_i >= kseqs_size) {
          kseqs_size *= 2;
          kseqs = (kseq_t ** )realloc(kseqs, kseqs_size*sizeof(kseq_t *));
        }
        kseqs[kseqs_i] = tmp;
        kseqs_i ++;
        size += tmp->seq.l;
      }
      if (kseqs_i >= g_mnmpOpt->mini_batch_size) {
        if (l_fragMode && kseqs[kseqs_i-1]->seq.l < CHECK_PAIR_THRES) {
          tmp_kseq = kseq_init(l_fp[0]->fp);
          while (kseq_read(tmp_kseq) >= 0) {
            if (mm_qname_same(tmp_kseq->name.s, kseqs[kseqs_i-1]->name.s)) {
              if (kseqs_i >= kseqs_size) {
                kseqs_size *= 2;
                kseqs = (kseq_t ** )realloc(kseqs, kseqs_size*sizeof(kseq_t *));
              }
              kseqs[kseqs_i] = tmp_kseq;
              kseqs_i ++;
              size += tmp_kseq->seq.l;
            }
            else {
              break;
            }
          }
        }
        break;
      }
    }

    if (kseqs_i == 0) {
      DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished KseqsRead"; 
      DLOG_IF(INFO, VLOG_IS_ON(2)) << "The last KseqsRead finished";
      break;
    }

    o_kseqsBatch.m_startSeqIdx = l_numProcessed;
    o_kseqsBatch.m_batchIdx    = l_batchCounter++;
    o_kseqsBatch.m_numSeq      = kseqs_i;
    o_kseqsBatch.kseqs         = kseqs;
    l_numProcessed += kseqs_i;
    
    DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished KseqsRead";
    DLOG_IF(INFO, VLOG_IS_ON(2)) << "number of seqs in o_kseqsBatch " << o_kseqsBatch.m_numSeq;
    DLOG_IF(INFO, VLOG_IS_ON(2)) << "m_startSeqIdx of o_kseqsBatch " << o_kseqsBatch.m_startSeqIdx;
    DLOG_IF(INFO, VLOG_IS_ON(2)) << "m_batchIdx of o_kseqsBatch " << o_kseqsBatch.m_batchIdx;
    
    this->pushOutput(o_kseqsBatch);
  }

  for (int f_i = 0; f_i < m_numFp; f_i++) {
    mm_bseq_close(l_fp[f_i]);
  }
  free(l_fp);

  return;
}