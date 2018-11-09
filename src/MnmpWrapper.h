#ifndef MNMP_FLOW_MINIMAP2_WRAPPER_H
#define MNMP_FLOW_MINIMAP2_WRAPPER_H

#include <string>

#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"

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
}

/* reimplement */
/* format.c */
std::string fc_write_sam_hdr(const mm_idx_t *idx, const std::string rg, const std::string ver, const std::string cmd);
/* map.c */
void fc_map_frag(const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, int *n_regs, mm_reg1_t **regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname);

#endif
