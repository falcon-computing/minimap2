#include "SeqsWriteStage.h"

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


int SeqsWrite::compute(BamsBatch const &i_bamsBatch) {
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started SeqsWrite";
  static char *l_writeMode[] = { "wb", "wb0", "w" };

  // Output File
  std::stringstream l_outputFname;
  samFile          *l_bamOutputFp;

  l_outputFname << FLAGS_output_dir << "/part-"
                << std::setw(6) << std::setfill('0') << i_bamsBatch.m_batchIdx << ".bam";
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
