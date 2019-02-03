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

const int max_kseq_buf_size = 100000; // set larger to prevent error when mini_batch_size is large

// copied from bseq.c
struct mm_bseq_file_s {
	gzFile fp;
	kseq_t *ks;
	mm_bseq1_t s;
};

// move kseq
void kseq_mov(kseq_t* ks1, kseq_t * ks2) {
  ks2->name = ks1->name;
  ks2->comment = ks1->comment;
  ks2->seq = ks1->seq;
  ks2->qual = ks1->qual;
  ks2->last_char = ks1->last_char;
  ks2->f = ks1->f;
  free(ks1);
}

KseqsRead::KseqsRead(int i_numFp, char **i_fn):kestrelFlow::SourceStage<KseqsBatch, INPUT_DEPTH>(),
                                              m_numFp(i_numFp),
                                              m_fn(i_fn)  
{
  for (int i = 0; i < INPUT_DEPTH; i++) {
    kseq_t * ks_new = (kseq_t*)calloc(max_kseq_buf_size, sizeof(kseq_t));
    kseq_buf buf_new;
    buf_new.ks = ks_new;
    buf_new.size = max_kseq_buf_size;
    kseq_queue.push(buf_new);
  }
}

KseqsRead::~KseqsRead() {
  for (int i = 0; i < INPUT_DEPTH; i++) {
    if (kseq_queue.empty()) {
      DLOG(ERROR) << "Internal memory corruption. No Kseqs_buf allocated.";
      break;
    }
    kseq_buf buf_tmp;
    kseq_queue.pop(buf_tmp);
    for (i = 0; i < buf_tmp.size; i++) {
      if (buf_tmp.ks[i].name.l) {
        free(buf_tmp.ks[i].name.s);
      }
      if (buf_tmp.ks[i].comment.l) {
        free(buf_tmp.ks[i].comment.s);
      }
      if (buf_tmp.ks[i].seq.l) {
        free(buf_tmp.ks[i].seq.s);
      }
      if (buf_tmp.ks[i].qual.l) {
        free(buf_tmp.ks[i].qual.s);
      }
    }
    free(buf_tmp.ks);
  }
}

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

  // kstream_t for all kseqs
  kstream_t ** kss = (kstream_t**)malloc(m_numFp*sizeof(kstream_t*));
  for (int i = 0; i < m_numFp; i++) {
    kss[i] = ks_init(l_fp[i]->fp);
  }

  while (true) {
    DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started KseqsRead";
    KseqsBatch o_kseqsBatch;

    kseq_buf buf_tmp;
    kseq_queue.pop(buf_tmp);

    kseq_t * kseqs = buf_tmp.ks;
    if (kseqs == NULL) {
      throw("kseq_buf is empty, shouldn't be.\n");
    }
    int kseqs_i = 0;
    int64_t size = 0;

    if (m_numFp > 1) {
      // read as in function mm_bseq_read_frag2
      while(1) {
        int n_read = 0;
        for (int i = 0; i < m_numFp; i++) {
          if (kseqs_i >= buf_tmp.size) {
              buf_tmp.size *= 2;
              kseqs = (kseq_t *)realloc(kseqs, buf_tmp.size * sizeof(kseq_t));
          }
          kseqs[kseqs_i].f = kss[i];
          if (kseq_read(kseqs+kseqs_i) >= 0) {
            n_read ++;
            kseqs_i ++;
            size += kseqs[kseqs_i].seq.l;
          }
        }
        if (n_read < m_numFp) {
          if (n_read > 0) {
              DLOG_IF(INFO, VLOG_IS_ON(2)) << "query files have different number of records; extra records skipped.";
          }
          break;
        }
        DLOG_IF(INFO, VLOG_IS_ON(4)) << "size " << size << " g_mnmpOpt->mini_batch_size " << g_mnmpOpt->mini_batch_size;
        if (size >= g_mnmpOpt->mini_batch_size) {
          break;
        }
      }
    }
    else {
      // read as in function mm_bseq_read3
      if (tmp_kseq != NULL) {
        kseq_mov(tmp_kseq, kseqs); // should always has 1 space in the beginning
        size += kseqs[0].seq.l;
        tmp_kseq = NULL;
      }
      while (1) {
        if (kseqs_i >= buf_tmp.size) {
          buf_tmp.size *= 2;
          kseqs = (kseq_t * )realloc(kseqs, buf_tmp.size*sizeof(kseq_t));
        }
        kseqs[kseqs_i].f = kss[0];
        if (kseq_read(kseqs + kseqs_i) < 0) {
          break;
        }
        kseqs_i ++;
        size += kseqs[kseqs_i].seq.l;
      }
      if (kseqs_i >= g_mnmpOpt->mini_batch_size) {
        if (l_fragMode && kseqs[kseqs_i-1].seq.l < CHECK_PAIR_THRES) {
          tmp_kseq = kseq_init(l_fp[0]->fp);
          while (kseq_read(tmp_kseq) >= 0) {
            if (mm_qname_same(tmp_kseq->name.s, kseqs[kseqs_i-1].name.s)) {
              if (kseqs_i >= buf_tmp.size) {
                buf_tmp.size *= 2;
                kseqs = (kseq_t * )realloc(kseqs, buf_tmp.size*sizeof(kseq_t));
              }
              kseq_mov(tmp_kseq, kseqs+kseqs_i);
              kseqs_i ++;
              size += kseqs[kseqs_i].seq.l;
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
    o_kseqsBatch.kseqs_buf     = buf_tmp;
    l_numProcessed += kseqs_i;
    
    DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished KseqsRead";
    DLOG_IF(INFO, VLOG_IS_ON(2)) << "number of seqs in o_kseqsBatch " << o_kseqsBatch.m_numSeq;
    DLOG_IF(INFO, VLOG_IS_ON(2)) << "m_startSeqIdx of o_kseqsBatch " << o_kseqsBatch.m_startSeqIdx;
    DLOG_IF(INFO, VLOG_IS_ON(2)) << "m_batchIdx of o_kseqsBatch " << o_kseqsBatch.m_batchIdx;
    
    this->pushOutput(o_kseqsBatch);
  }

  for (int f_i = 0; f_i < m_numFp; f_i++) {
    ks_destroy(kss[f_i]);
    mm_bseq_close(l_fp[f_i]);
  }
  free(l_fp);

  return;
}