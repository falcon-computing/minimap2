#include "MarkDupStage2.h"
#include <cstring>

#define DUP_FLAG 1024
#define NP_FLAG 256

static splitLine_t * bamToSplitLine(bam_hdr_t* hdr, bam1_t* bam_record) {
  splitLine* sline = getSplitLine();
  kstring_t ks = { 0, 0, NULL };
  sam_format1(hdr, bam_record, &ks);
  sline->bufLen = ks.l;
  strcpy(sline->buffer, ks.s);
  free(ks.s);
  splitSplitLine(sline, 12);
  return sline;
}

static bool checkSplitLineDup(splitLine_t * sline) {
  return (bool)(sline->flag & DUP_FLAG);
}

static void markDupBam(bam1_t * bam) {
  bam->core.flag = (bam->core.flag | DUP_FLAG);
  return;
}

void MarkDupStage::InitializeState(bam_hdr_t* hdr) {
  state = makeState();
  state->seqLens = (UINT32*)calloc(1, sizeof(UINT32));
  state->seqOffs = (UINT64*)calloc(1, sizeof(UINT64));
  state->seqs[strdup("*")] = 0;
  state->seqLens[0] = padLength(0);
  state->seqOffs[0] = 0;
#ifdef USE_HTSLIB
  UINT64 totalLen = 0;
  for(int i = 0; i < hdr->n_targets; i++) {
    char * seqID = hdr->target_name[i];;
    UINT32 seqLen = (UINT32)hdr->target_len[i];
    UINT64 seqOff = totalLen;
    totalLen += (UINT64)(seqLen + 1);
    if(i % 32768 == 1) {
      state->seqLens = (UINT32*)realloc(state->seqLens, (i + 32768)*sizeof(UINT32));
      state->seqOffs = (UINT64*)realloc(state->seqOffs, (i + 32768)*sizeof(UINT64));
    }
    state->seqs[strdup(seqID)] = i;
    state->seqLens[i] = seqLen;
    state->seqOffs[i] = seqOff;
  }
  int binCount = (totalLen >> BIN_SHIFT);
  if (binCount >= (1 << 15)) {
    //Error Too many sequences in header of input sam file.
  }
  state->binCount = binCount;
  state->sigArraySize = (binCount * 2 + 1) * (binCount * 2 + 1) + 1;
  //state->sigs = (sigSet_t *) malloc(state->sigArraySize * sizeof(sigSet_t));
  state->sigs = (sigSet_t *) new sigSet_t[state->sigArraySize];
  if (state->sigs == NULL) fatalError("samblaster: Unable to allocate signature set array.");
  //for (UINT32 i=0; i<state->sigArraySize; i++) state->sigs[i].hashTableInit(); 
#else
#endif
}

BamsBatch MarkDupStage::compute(BamsBatch const & input) {
    DLOG(INFO) << "Started MarkDup()";
    int batch_num = input.m_numBams;
    splitLine_t* line = bamToSplitLine(hdr_, input.m_bams[0]);
    splitLine_t* head = line;
    int count = 1;
    splitLine_t* last = line;
    for (int i = 1; i < batch_num; i++) {
      splitLine_t* nextLine =  bamToSplitLine(hdr_, input.m_bams[i]);
      if (nextLine == NULL) {
        break;
      }
      if (strcmp(line->fields[QNAME], nextLine->fields[QNAME]) != 0) {
        mtx_.lock();
        markDupsDiscordants(line, state);
        mtx_.unlock();
        splitLine_t* tmp = line;
        for (int j = 0; j < count; j++) {
          if (checkSplitLineDup(tmp)) {
            markDupBam(input.m_bams[i- count +j]);
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
    markDupsDiscordants(line, state);
    mtx_.unlock();
    splitLine_t* tmp = line;
    for(int j = 0; j < count; j++) {
      if (checkSplitLineDup(tmp)) {
        markDupBam(input.m_bams[batch_num - count + j]);
      }
      splitLine* release = tmp;
      tmp = tmp->next;
      deleteSplitLine(release);
    }
    DLOG(INFO) << "Finished MarkDup()";
    return input;
}
