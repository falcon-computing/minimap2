#ifndef MNMP_FLOW_MINIMAP2_WRAPPER_H
#define MNMP_FLOW_MINIMAP2_WRAPPER_H

#include <string>

#include "bseq.h"
#include "ksw2.h"
#include "minimap.h"
#include "mmpriv.h"

#include "MnmpGlobal.h"
#include "MnmpOptions.h"
#include "MnmpData.h"
/* static export */
extern "C" {
  /* map.c */
  mm_bseq_file_t **open_bseqs(int n, const char **fn);
  void collect_minimizers(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, mm128_v *mv);
  mm128_t *collect_seed_hits_heap(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi, const char *qname, const mm128_v *mv, int qlen, int64_t *n_a, int *rep_len,
                                  int *n_mini_pos, uint64_t **mini_pos);
  mm128_t *collect_seed_hits(void *km, const mm_mapopt_t *opt, int max_occ, const mm_idx_t *mi, const char *qname, const mm128_v *mv, int qlen, int64_t *n_a, int *rep_len,
                             int *n_mini_pos, uint64_t **mini_pos);
  void chain_post(const mm_mapopt_t *opt, int max_chain_gap_ref, const mm_idx_t *mi, void *km, int qlen, int n_segs, const int *qlens, int *n_regs, mm_reg1_t *regs, mm128_t *a);
  mm_reg1_t *align_regs(const mm_mapopt_t *opt, const mm_idx_t *mi, void *km, int qlen, const char *seq, int *n_regs, mm_reg1_t *regs, mm128_t *a);

  /* align.c */
  mm_reg1_t *mm_insert_reg(const mm_reg1_t *r, int i, int *n_regs, mm_reg1_t *regs);
  void mm_align1(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, uint8_t *qseq0[2], mm_reg1_t *r, mm_reg1_t *r2, int n_a, mm128_t *a, ksw_extz_t *ez, int splice_flag);
  int mm_align1_inv(void *km,const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, uint8_t *qseq0[2], const mm_reg1_t *r1, const mm_reg1_t *r2, mm_reg1_t *r_inv, ksw_extz_t *ez);
}

/* reimplement */
/* format.c */
std::string fc_write_sam_hdr(const mm_idx_t *idx, const std::string rg, const std::string ver, const std::string cmd);
/* map.c */
void fc_map_frag_chain(const mm_mapopt_t *opt, const mm_idx_t *mi, const int i_fragIdx, int n_segs, const int *qlens, const char **seqs, int *n_regs, mm_reg1_t **regs, const char *qname, fragExtSOA *io_fragExtSOA, int *o_repLen, int *o_fragGap);
void fc_map_frag_align(const mm_mapopt_t *opt, const mm_idx_t *mi, const int i_fragIdx, int n_segs, const int *qlens, const char **seqs, int *n_regs, mm_reg1_t **regs, const char *qname, const int i_repLen, const int i_fragGap, fragExtSOA *io_fragExtSOA);
mm_tbuf_t *fc_tbuf_init(void);
void fc_tbuf_destroy(mm_tbuf_t *);

#endif
