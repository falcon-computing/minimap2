#include "MarkDupStage.h"
#include <cstring>

#include "MnmpGlobal.h"
#include "MnmpOptions.h"
#include "MnmpWrapper.h"
#include "MnmpUtils.h"
#include "minimap.h"
#include "mmpriv.h"
#include "kalloc.h"

#include "kflow/Common.h"
#include "kflow/Pipeline.h"
#include "kflow/SourceStage.h"
#include "kflow/SinkStage.h"
#include "samblaster.h"
#include "sbhash.h"

#define DUP_FLAG 1024
#define NP_FLAG 256

static bool bam1_lt(const bam1_t *a, const bam1_t *b) {
  return ((uint64_t)a->core.tid<<32|(a->core.pos+1)<<1|bam_is_rev(a))
       < ((uint64_t)b->core.tid<<32|(b->core.pos+1)<<1|bam_is_rev(b));
}

static splitLine_t * bamToSplitLine(bam_hdr_t* head, bam1_t* bam_record) {
  splitLine* sline = getSplitLine();
  kstring_t ks = { 0, 0, NULL };
  sam_format1(head, bam_record, &ks);
  sline->bufLen = ks.l;
  strcpy(sline->buffer, ks.s);
  free(ks.s);
  splitSplitLine(sline, 12);
  return sline;
}

static splitLine_t * kstringToSplitLine(kstring_t * ks) {
  splitLine* sline = getSplitLine();
  sline->bufLen = ks->l;
  memcpy(sline->buffer, ks->s, ks->l*sizeof(char));

  splitSplitLine(sline, 12);
  return sline;
}

// static splitLine_t * readSeq(bam_hdr_t* head, bseq1_t seq) {
//   if (seq.bams->bams[0] != NULL) {
//     if (seq.bams->bams[0]->core.flag & NP_FLAG) {
//       LOG(INFO)<<"alignment not primary.";
//     }
//     splitLine* line = bamToSplitLine(head, seq.bams->bams[0]);
//     if (line->bufLen < 1) {
//       return NULL;
//     }
//     return line;
//   }
//   return NULL;
// }

static bool checkSplitLineDup(splitLine_t * sline) {
  return (bool)(sline->flag & DUP_FLAG);
}

// static void markDupSeq(bseq1_t* seq) {
//   for(int i = 0; i < seq->bams->l; i++) {
//     if (i >= 2) {
//       DLOG(INFO)<<"more than 1 alignments marked: "<<i;
//     }
//     if(seq->bams->bams[i] != NULL) {
//       seq->bams->bams[i]->core.flag = (seq->bams->bams[i]->core.flag | DUP_FLAG);
//     }
//   }
//   return;
// }

static void markDupBam1t(bam1_t* bam) {
  bam->core.flag = (bam->core.flag | DUP_FLAG);
  return;
}

void MarkDupStage::InitializeState(bam_hdr_t* head) {
  state_ = makeState();
  state_->seqLens = (UINT32*)calloc(1, sizeof(UINT32));
  state_->seqOffs = (UINT64*)calloc(1, sizeof(UINT64));
  state_->seqs[strdup("*")] = 0;
  state_->seqLens[0] = padLength(0);
  state_->seqOffs[0] = 0;
#ifdef USE_HTSLIB
  UINT64 totalLen = 0;
  for(int i = 0; i < head->n_targets; i++) {
    char * seqID = head->target_name[i];;
    UINT32 seqLen = (UINT32)head->target_len[i];
    UINT64 seqOff = totalLen;
    totalLen += (UINT64)(seqLen + 1);
    if(i % 32768 == 1) {
      state_->seqLens = (UINT32*)realloc(state_->seqLens, (i + 32768)*sizeof(UINT32));
      state_->seqOffs = (UINT64*)realloc(state_->seqOffs, (i + 32768)*sizeof(UINT64));
    }
    state_->seqs[strdup(seqID)] = i;
    state_->seqLens[i] = seqLen;
    state_->seqOffs[i] = seqOff;
  }
  int binCount = (totalLen >> BIN_SHIFT);
  if (binCount >= (1 << 15)) {
    //Error Too many sequences in header of input sam file.
  }
  state_->binCount = binCount;
  state_->sigArraySize = (binCount * 2 + 1) * (binCount * 2 + 1) + 1;
  //state_->sigs = (sigSet_t *) malloc(state_->sigArraySize * sizeof(sigSet_t));
  state_->sigs = (sigSet_t *) new sigSet_t[state_->sigArraySize];
  if (state_->sigs == NULL) fatalError("samblaster: Unable to allocate signature set array.");
  //for (UINT32 i=0; i<state_->sigArraySize; i++) state_->sigs[i].hashTableInit(); 
#else
#endif
}

BamsBatch MarkDupStage::compute(AlignsBundle const & input) {
  bam_hdr_t *l_bamHeader = g_bamHeader;
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started Markdup";

  int l_numBamsEst = 0;
  for (int l_ba = 0; l_ba < input.m_batches->size(); l_ba++) {
    AlignsBatch &l_alignsBatch = (*input.m_batches)[l_ba];
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
  for (int l_ba = 0; l_ba < input.m_batches->size(); l_ba++) {
    AlignsBatch &l_alignsBatch = (*input.m_batches)[l_ba];
    for (int l_fr = 0; l_fr < l_alignsBatch.m_numFrag; l_fr++) {
      int l_segOffset = l_alignsBatch.m_segOff[l_fr];
      int l_numSegs   = l_alignsBatch.m_numSeg[l_fr];
      
      splitLine_t * splitLines = NULL;
      splitLine_t * next_splitLine = NULL;
      std::vector<bam1_t*> l_bams;
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
            mm_write_sam2(&l_samStrBuf, g_minimizer, l_seq, l_sg-l_segOffset, l_rg, 
                          l_numSegs, l_numRegsArr, l_regsArr, NULL, g_mnmpOpt->flag);
            if (l_numBams >= l_numBamsEst) {
              l_numBamsEst += (int)(0.2*l_numBamsEst);
              l_bamsArr = (bam1_t**)realloc(l_bamsArr, l_numBamsEst*sizeof(bam1_t*));
            }
            bam1_t *l_bam = bam_init1();
            l_retCode = sam_parse1(&l_samStrBuf, l_bamHeader, l_bam);
            next_splitLine = bamToSplitLine(l_bamHeader, l_bam);
            if (splitLines == NULL) {
              splitLines = next_splitLine;
            }
            else {
              splitLine_t * tmp = splitLines;
              while (tmp->next != NULL) {
                tmp = tmp->next;
              }
              tmp->next = next_splitLine;
            }
            
            l_bams.push_back(l_bam);
            l_bamsArr[l_numBams++] = l_bam;
          }
        }
        else if (g_mnmpOpt->flag & (MM_F_OUT_SAM|MM_F_PAF_NO_HIT)) { // output an empty hit, if requested
          mm_write_sam2(&l_samStrBuf, g_minimizer, l_seq, l_sg-l_segOffset, -1, 
                        l_numSegs, l_numRegsArr, l_regsArr, NULL, g_mnmpOpt->flag);
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
      if (splitLines != NULL) {
        mtx_.lock();
        markDupsDiscordants(splitLines, state_);
        mtx_.unlock();
        if (checkSplitLineDup(splitLines)) {
          for (int ii = 0; ii < l_bams.size(); ii++) {
            markDupBam1t(l_bams[ii]);
          }
        }
      }
      splitLine_t * tmp = splitLines;
      while (tmp != NULL) {
        splitLine_t * tmpp = tmp->next;
        deleteSplitLine(tmp);
        tmp = tmpp;
      }
    } // close loop over fragments
    // Free a batch of alignment
    free(l_alignsBatch.m_reg);
    free(l_alignsBatch.m_numReg);
    free(l_alignsBatch.m_seqs); // seg_off, n_seg, rep_len and frag_gap were allocated with reg; no memory leak here
  } // close loop over batches
  delete input.m_batches;
  l_bamsArr = (bam1_t**)realloc(l_bamsArr, l_numBams*sizeof(bam1_t*));
  free(l_samStrBuf.s);
  if (!FLAGS_disable_sort && FLAGS_disable_bucketsort) {
    std::sort(l_bamsArr, l_bamsArr+l_numBams, bam1_lt);
  }
  BamsBatch o_bamsBatch;
  o_bamsBatch.m_batchIdx = input.m_bundleIdx;
  o_bamsBatch.m_bams     = l_bamsArr;
  o_bamsBatch.m_numBams  = l_numBams;
  
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished Markdup";

  return o_bamsBatch;
#if 0
//uint64_t read_seq_time = 0;
//uint64_t mark_dup_time = 0;
    //uint64_t all_start = getUs();
    DLOG(INFO) << "Started MarkDup()";
    //uint64_t read_seq_s = getUs();
    //for SeqsReord input
    std::vector<SeqsRecord> tmp_vec;
    tmp_vec.push_back(input);
    std::vector<SeqsRecord> * seqsRecord = &tmp_vec;
    //uint64_t read_seq_e = getUs();
    //read_seq_time += (read_seq_e - read_seq_s);
    for (int k =0; k < (*seqsRecord).size(); k++) {
      int batch_num = (*seqsRecord)[k].batch_num;
      splitLine_t* line = readSeq(head, (*seqsRecord)[k].seqs[0]);
      if (line != NULL) {
        //DLOG(INFO) << "fisrt read is not NULL.";
      }
      splitLine_t* head = line;
      int count = 1;
      splitLine_t* last = line;
      for (int i = 1; i < batch_num; i++) {
        //read_seq_s = getUs();
        splitLine_t* nextLine =  readSeq(head, (*seqsRecord)[k].seqs[i]);
        //read_seq_e = getUs();
        //read_seq_time += (read_seq_e - read_seq_s);
        if (nextLine == NULL) {
          break;
        }
        if (strcmp(line->fields[QNAME], nextLine->fields[QNAME]) != 0) {
          //DLOG(INFO) << "before md " << line->fields[QNAME] << " " << nextLine->fields[QNAME];
          mtx_.lock();
          markDupsDiscordants(line, state_);
          mtx_.unlock();
          splitLine_t* tmp = line;
          for (int j = 0; j < count; j++) {
            if (checkSplitLineDup(tmp)) {
              //uint64_t mark_dup_s = getUs();
              markDupSeq(&((*seqsRecord)[k].seqs[i- count +j]));
              //uint64_t mark_dup_e = getUs();
              //mark_dup_time += (mark_dup_e - mark_dup_s);
            }
            splitLine* release = tmp;
            tmp = tmp->next;
            deleteSplitLine(release);
          }
          last = line = nextLine;
          count = 1;
        }
        else{
          last->next = nextLine;
          last = nextLine;
          count += 1;
        }
      } 

      mtx_.lock();
      markDupsDiscordants(line, state_);
      mtx_.unlock();
      splitLine_t* tmp = line;
      for(int j = 0; j < count; j++) {
        if (checkSplitLineDup(tmp)) {
          //uint64_t mark_dup_s = getUs();
          markDupSeq(&((*seqsRecord)[k].seqs[batch_num - count + j]));
          //uint64_t mark_dup_e = getUs();
          //mark_dup_time += (mark_dup_e - mark_dup_s);
        }
        splitLine* release = tmp;
        tmp = tmp->next;
        deleteSplitLine(release);
      }
    }
    //pushOutput(input); // for mapPartitionStage
    DLOG(INFO) << "Finished MarkDup()";
    //uint64_t all_end = getUs();
    //uint64_t all_diff = all_end - all_start;
    //DLOG(INFO) << "MdStage AllTime = " << all_diff;
    //DLOG(INFO) << "ReadSeq Time = " << read_seq_time;
    //DLOG(INFO) << "MarkDup Time = " << mark_dup_time;
    return input;
  // } // for MapPartitionStage
#if 0
  tmp = head;
  while(tmp->next != NULL) {
    splitLine_t* release = tmp;
    tmp = tmp->next;
    deleteSplitLine(release);
  }
  deleteSplitLine(tmp);
#endif
  //mtx_.unlock();
#endif
}
