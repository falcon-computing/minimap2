#include "MnmpWrapper.h"

#include <string>
#include <fstream>
#include <sstream>

#include <glog/logging.h>

#include "kalloc.h"
#include "khash.h"
#include "ksw2.h"

std::string fc_write_sam_hdr(const mm_idx_t *idx, const std::string rg, const std::string ver, const std::string cmd)
{
  //kstring_t str = {0,0,0};
  std::stringstream l_hdrStr;
  if (idx) {
    for (uint32_t i = 0; i < idx->n_seq; ++i)
      l_hdrStr << "@SQ\tSN:" << idx->seq[i].name << "\tLN:" << idx->seq[i].len << "\n";
      //mm_sprintf_lite(&str, "@SQ\tSN:%s\tLN:%d\n", idx->seq[i].name, idx->seq[i].len);
  }
  if (rg.length() > 0) {
    l_hdrStr << rg;
    //sam_write_rg_line(&str, rg);
  }
  l_hdrStr << "@PG\tID:minimap2\tPN:minimap2";
  //mm_sprintf_lite(&str, "@PG\tID:minimap2\tPN:minimap2");
  if (ver.length() > 0) {
    l_hdrStr << "\tVN:" << ver;
    //mm_sprintf_lite(&str, "\tVN:%s", ver);
  }
  if (cmd.length() > 0) {
    l_hdrStr << "\tCL:" << cmd << std::endl;
    //int i;
    //mm_sprintf_lite(&str, "\tCL:minimap2");
    //for (i = 1; i < argc; ++i)
    //  mm_sprintf_lite(&str, " %s", argv[i]);
  }
  return l_hdrStr.str();
  //mm_err_puts(str.s);
  //free(str.s);
}

mm_tbuf_t *fc_tbuf_init(void)
{
  mm_tbuf_t *b;
  b = (mm_tbuf_t*)calloc(1, sizeof(mm_tbuf_t));
  b->km = NULL;
  return b;
}

void fc_tbuf_destroy(mm_tbuf_t *b)
{
  if (b == NULL) return;
  free(b);
}

void fc_map_frag_chain(const mm_mapopt_t *opt, const mm_idx_t *mi, const int i_fragIdx, int n_segs, const int *qlens, const char **seqs, int *n_regs, mm_reg1_t **regs, const char *qname, fragExtSOA *io_fragExtSOA, int *o_repLen, int *o_fragGap)
{
	int i, j, n_mini_pos;
	int max_chain_gap_qry, max_chain_gap_ref, is_splice = !!(opt->flag & MM_F_SPLICE), is_sr = !!(opt->flag & MM_F_SR);
	int64_t n_a;

        uint32_t l_hash;
        mm_reg1_t *l_regs0Arr;
        int l_numRegs0;
        uint64_t *l_uArr;
        uint64_t *l_miniPosArr;
        mm128_t *l_anchorArr;
        mm128_v l_mv = {0, 0, 0};

        int qlen_sum = 0;
	for (int l_sg = 0; l_sg < n_segs; ++l_sg)
		qlen_sum += qlens[l_sg], n_regs[l_sg] = 0, regs[l_sg] = 0;

	if (qlen_sum == 0 || n_segs <= 0 || n_segs > MM_MAX_SEG) return;

	l_hash  = qname? __ac_X31_hash_string(qname) : 0;
	l_hash ^= __ac_Wang_hash(qlen_sum) + __ac_Wang_hash(opt->seed);
	l_hash  = __ac_Wang_hash(l_hash);

	collect_minimizers(NULL, opt, mi, n_segs, qlens, seqs, &l_mv);
	if (opt->flag & MM_F_HEAP_SORT) l_anchorArr = collect_seed_hits_heap(NULL, opt, opt->mid_occ, mi, qname, &l_mv, qlen_sum, &n_a, o_repLen, &n_mini_pos, &l_miniPosArr);
	else l_anchorArr = collect_seed_hits(NULL, opt, opt->mid_occ, mi, qname, &l_mv, qlen_sum, &n_a, o_repLen, &n_mini_pos, &l_miniPosArr);

#if 0
        for (int l_a = 0; l_a < n_a; l_a++)
          DLOG(INFO) << "SD\t" << mi->seq[l_anchorArr[l_a].x<<1>>33].name << "\t"
                     << (int32_t)l_anchorArr[l_a].x << "\t" << "+-"[l_anchorArr[l_a].x>>63] << "\t"
                     << (int32_t)l_anchorArr[l_a].y << "\t" << (int32_t)(l_anchorArr[l_a].y>>32&0xff) << "\t"
                     << (l_a == 0) ? 0: ((int32_t)l_anchorArr[l_a].y-(int32_t)l_anchorArr[l_a-1].y) - ((int32_t)l_anchorArr[l_a].x-(int32_t)l_anchorArr[l_a-1].x);
#endif

	// set max chaining gap on the query and the reference sequence
	if (is_sr)
		max_chain_gap_qry = qlen_sum > opt->max_gap? qlen_sum : opt->max_gap;
	else max_chain_gap_qry = opt->max_gap;
	if (opt->max_gap_ref > 0) {
		max_chain_gap_ref = opt->max_gap_ref; // always honor mm_mapopt_t::max_gap_ref if set
	} else if (opt->max_frag_len > 0) {
		max_chain_gap_ref = opt->max_frag_len - qlen_sum;
		if (max_chain_gap_ref < opt->max_gap) max_chain_gap_ref = opt->max_gap;
	} else max_chain_gap_ref = opt->max_gap;

	l_anchorArr = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, l_anchorArr, &l_numRegs0, &l_uArr, NULL);

	if (opt->max_occ > opt->mid_occ && *o_repLen > 0) {
		int rechain = 0;
		if (l_numRegs0 > 0) { // test if the best chain has all the segments
			int n_chained_segs = 1, max = 0, max_i = -1, max_off = -1, off = 0;
			for (i = 0; i < l_numRegs0; ++i) { // find the best chain
				if (max < (int)(l_uArr[i]>>32)) max = l_uArr[i]>>32, max_i = i, max_off = off;
				off += (uint32_t)l_uArr[i];
			}
			for (i = 1; i < (int32_t)l_uArr[max_i]; ++i) // count the number of segments in the best chain
				if ((l_anchorArr[max_off+i].y&MM_SEED_SEG_MASK) != (l_anchorArr[max_off+i-1].y&MM_SEED_SEG_MASK))
					++n_chained_segs;
			if (n_chained_segs < n_segs)
				rechain = 1;
		} else rechain = 1;
		if (rechain) { // redo chaining with a higher max_occ threshold
			kfree(NULL, l_anchorArr);
			kfree(NULL, l_uArr);
			kfree(NULL, l_miniPosArr);
			if (opt->flag & MM_F_HEAP_SORT) l_anchorArr = collect_seed_hits_heap(NULL, opt, opt->max_occ, mi, qname, &l_mv, qlen_sum, &n_a, o_repLen, &n_mini_pos, &l_miniPosArr);
			else l_anchorArr = collect_seed_hits(NULL, opt, opt->max_occ, mi, qname, &l_mv, qlen_sum, &n_a, o_repLen, &n_mini_pos, &l_miniPosArr);
			l_anchorArr = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, l_anchorArr, &l_numRegs0, &l_uArr, NULL);
		}
	}

	l_regs0Arr = mm_gen_regs(NULL, l_hash, qlen_sum, l_numRegs0, l_uArr, l_anchorArr);

#if 0
        for (int l_rg0; l_rg0 < l_numRegs0; ++l_rg0)
          for (int l_cn = l_regs0Arr[l_rg0].as; l_cn < l_regs0Arr[l_rg0].as + l_regs0Arr[l_rg0].cnt; ++l_cn)
            DLOG(INFO) << "CN\t" << l_rg0 << "\t" << mi->seq[l_anchorArr[l_cn].x<<1>>33].name << "\t"
                       << (int32_t)l_anchorArr[l_cn].x << "\t" << "+-"[l_anchorArr[l_cn].x>>63] << "\t"
                       << (int32_t)l_anchorArr[l_cn].y << "\t" << (int32_t)(l_anchorArr[l_cn].y>>32&0xff) << "\t"
                       << (l_cn == l_regs0Arr[l_rg0].as) ? 0
                                  : ((int32_t)l_anchorArr[l_cn].y-(int32_t)l_anchorArr[l_cn-1].y) - ((int32_t)l_anchorArr[l_cn].x-(int32_t)l_anchorArr[l_cn-1].x);
#endif

	chain_post(opt, max_chain_gap_ref, mi, NULL, qlen_sum, n_segs, qlens, &l_numRegs0, l_regs0Arr, l_anchorArr);
	if (!is_sr) mm_est_err(mi, qlen_sum, l_numRegs0, l_regs0Arr, l_anchorArr, n_mini_pos, l_miniPosArr);

        *o_fragGap = max_chain_gap_ref;
        //*o_repLen  = rep_len;

        io_fragExtSOA->m_regs0Arr[i_fragIdx] = l_regs0Arr;
        io_fragExtSOA->m_numRegs0[i_fragIdx] = l_numRegs0;

        io_fragExtSOA->m_hash[i_fragIdx] = l_hash;

        io_fragExtSOA->m_anchorArr[i_fragIdx] = l_anchorArr; // numAnchor is not used in map_align

        kfree(NULL, l_mv.a);
        kfree(NULL, l_uArr);
        kfree(NULL, l_miniPosArr);
}

void fc_map_frag_align(const mm_mapopt_t *opt, const mm_idx_t *mi, const int i_fragIdx, int n_segs, const int *qlens, const char **seqs, int *n_regs, mm_reg1_t **regs, const char *qname, const int i_repLen, const int i_fragGap, fragExtSOA *io_fragExtSOA)
{
        uint32_t l_hash = io_fragExtSOA->m_hash[i_fragIdx];
        mm_reg1_t *l_regs0Arr = io_fragExtSOA->m_regs0Arr[i_fragIdx];
        int l_numRegs0 = io_fragExtSOA->m_numRegs0[i_fragIdx];
        mm128_t *l_anchorArr = io_fragExtSOA->m_anchorArr[i_fragIdx];

	int is_sr = !!(opt->flag & MM_F_SR);

	if (n_segs == 1) { // uni-segment
		l_regs0Arr = align_regs(opt, mi, NULL, qlens[0], seqs[0], &l_numRegs0, l_regs0Arr, l_anchorArr);
		mm_set_mapq(NULL, l_numRegs0, l_regs0Arr, opt->min_chain_score, opt->a, i_repLen, is_sr);
		n_regs[0] = l_numRegs0, regs[0] = l_regs0Arr;
	} else { // multi-segment
		mm_seg_t *seg;
		seg = mm_seg_gen(NULL, l_hash, n_segs, qlens, l_numRegs0, l_regs0Arr, n_regs, regs, l_anchorArr); // split fragment chain to separate segment chains
		free(l_regs0Arr);
		for (int i = 0; i < n_segs; ++i) {
			mm_set_parent(NULL, opt->mask_level, n_regs[i], regs[i], opt->a * 2 + opt->b, opt->flag&MM_F_HARD_MLEVEL); // update mm_reg1_t::parent
			regs[i] = align_regs(opt, mi, NULL, qlens[i], seqs[i], &n_regs[i], regs[i], seg[i].a);
			mm_set_mapq(NULL, n_regs[i], regs[i], opt->min_chain_score, opt->a, i_repLen, is_sr);
		}
		mm_seg_free(NULL, n_segs, seg);
		if (n_segs == 2 && opt->pe_ori >= 0 && (opt->flag&MM_F_CIGAR))
			mm_pair(NULL, i_fragGap, opt->pe_bonus, opt->a * 2 + opt->b, opt->a, qlens, n_regs, regs); // pairing
	}

	kfree(NULL, l_anchorArr);
}

bool packInput(int const                 i_fragId,
               int const                 i_numSegs,
               int const                 i_segOff,
               int const                 i_repLen,
               int const                 i_fragGap,
               int const                 i_numRegs0,
               mm_reg1_t * const         i_regs0,
               mm128_t   * const         i_anks,
               uint32_t    const         i_hash,
               std::vector<mm_seg_t*>   &o_segChnsInChunkArr,
               std::vector<int>         &o_seqReference,
               std::vector<align_input> &o_fpgaChunkInput)
{
  const int l_fr = io_fragId;
  const int l_numSegs = i_chainsBatch.m_numSeg[l_fr];
  const int l_segOff  = i_chainsBatch.m_segOff[l_fr];

  fragExtSOA * const l_fragExtSOA = i_chainsBatch.m_fragExtSOA;

  int l_queryLens[MM_MAX_SEG];
  const char *l_querySeqs[MM_MAX_SEG];
  for (int l_sg = 0; l_sg < l_numSegs; l_sg++) {
    l_queryLens[l_sg] = i_chainsBatch.m_seqs[l_segOff + l_sg].l_seq;
    l_querySeqs[l_sg] = i_chainsBatch.m_seqs[l_segOff + l_sg].seq;
  }

  mm_seg_t *l_segChns = NULL;
  if (l_numSegs > 1) {
    // split fragment chain to separate segment chains (original scope: mm_map_frag/fc_map_frag_align)
    l_segChns = mm_seg_gen(NULL,
                          l_hash,
                          l_numSegs,
                          l_queryLens, 
                          l_numRegs0,
                          l_regs0Arr,
                          l_numRegsOfSgmArr,
                          l_regsOfSgmArr,
                          l_anks);
    free(l_regs0Arr);

    // update mm_reg1_t::parent (original scope: mm_map_frag/fc_map_frag_align)
    for (int l_sg = 0; l_sg < l_numSegs; l_sg++) {
      mm_set_parent(NULL,
                    g_mnmpOpt->mask_level,
                    l_numRegsOfSgmArr[l_sg],
                    l_regsOfSgmArr[l_sg],
                    g_mnmpOpt->a * 2 + g_mnmpOpt->b,
                    g_mnmpOpt->flag&MM_F_HARD_MLEVEL);
    }
  }
  o_segChnsInChunkArr.push_back(l_segChns);

  for (int l_sg = 0; l_sg < l_numSegs; l_sg++) {
    /* the number of  regions of a segment in the fragment, aka a sequence in the pair */
    int         l_numRegs = l_numRegsOfSgmArr[l_sg];
    /* pointer to regions of a segment in the fragment, aka a sequence in the pair */
    mm_reg1_t * l_regs    = l_regsOfSgmArr[l_sg];

    // alias anchors (original scope: align_regs)
    mm128_t *l_anchors = (l_numSegs > 1) ? l_segChns[l_sg].a : l_anks;
    // get squeezed number of anchors (original scope: mm_align_skeleton, after query seqs encoding)
    int32_t l_numAnchors = mm_squeeze_a(NULL, l_numRegs, l_regs, l_anchors);

    if (l_queryLens[l_sg] > MAX_SEQ_LENGTH || l_numAnchors > MAX_ANKOR_NUM) {
      return false;
    }
    
    // encode the query sequence (original scope: mm_align_skeleton, before regions squeezing)
    uint8_t l_encodedQseq[2][MAX_SEQ_LENGTH];
    int  const         l_qlen = l_queryLens[l_sg];
    char const * const l_qseq = l_querySeqs[l_sg];
    for (int l_bs = 0; l_bs < l_qlen; l_bs++) {
      l_encodedQseq[0][l_bs] = seq_nt4_table[(uint8_t)l_qseq[l_bs]];
      l_encodedQseq[1][l_queryLens[l_sg]-l_bs-1] = (l_encodedQseq[0][l_bs] < 4) ? 3-l_encodedQseq[0][l_bs] : 4;
    }

    // generate align kernel inputs for all regions
    for (int l_rg = 0; l_rg < l_numRegs; l_rg++) {
      align_input l_newAlignInput;
      l_newAlignInput.qlen = l_qlen;
      for (int l_bs = 0; l_bs < l_qlen; l_bs++) {
        l_newAlignInput.qseq[0][l_bs] = l_encodedQseq[0][l_bs];
        l_newAlignInput.qseq[1][l_bs] = l_encodedQseq[1][l_bs];
      }

      l_newAlignInput.n_a = l_numAnchors;
      memcpy(l_newAlignInput.a, l_anchors, l_numAnchors*sizeof(mm128_t));

      /* leave region[1] as blank */
      l_newAlignInput.region[0].orig = l_regs[l_rg];
      /* leave dp_score[1] as blank */
      l_newAlignInput.dp_score[0] = l_regs[l_rg].p->dp_score;

      l_newAlignInput.data_size = sizeof(align_input) - 16 * MAX_ANKOR_NUM + 16 * l_numAnchors;

      o_seqReference.push_back(l_fr);
      o_seqReference.push_back(l_sg);
      o_seqReference.push_back(l_rg);

      o_alignInput.push_back(l_newAlignInput);
    } // close loop over regions
  } // close loop over segments

  return true;
}

bool needRecalculate(int           i_numOutputs,
                     align_output *i_alignOutputs)
{
  for (int l_i = 0; l_i < i_numOutputs; l_i++) {
    if ( i_alignOutputs[l_i].region[1].cnt > 0                  ||
         ( l_i > 0 && i_alignOutputs[l_i].region[0].split_inv )   ) {
      return true;
    }
  }
  return false;
}

void unpackOutput(int i_fragId,
                  int i_sgmId,
                  int i_regId,
                  align_input        i_alignInput,
                  align_output       i_alignOutput,
                  mm_reg1_t*         i_regs,
                  int const          i_numRegs,
                  mm_seg_t*          i_segs,
                  int const          i_numSegs,
                  ChainsBatch const &io_chainsBatch)
{
  const int l_fr = i_fragId; 
  const int l_sg = i_sgmId; 
  const int l_rg = i_regId; 

  const int l_segOff = i_chainsBatch.m_segOff[l_fr];
  const int l_numReg = i_chainsBatch.m_numReg[l_segOff+l_sg];

  // post-process on regions (original scope: mm_align_skeleton)
  if (l_rg+1 == l_numReg) {
    mm_filter_regs(g_mnmpOpt,
                   i_alignInput.qlen,
                  &i_chainsBatch.m_numReg[l_segOff+l_sg], 
                  &i_chainsBatch.m_reg[l_segOff+l_sg]);
    mm_hit_sort(NULL, &i_chainsBatch.m_numReg[l_segOff+l_sg], i_chainsBatch.m_reg[l_segOff+l_sg]);
  }

  // don't choose primary mapping(s) (original scope: align_regs)
  if (!(opt->flag & MM_F_ALL_CHAINS)) {
    mm_set_parent(NULL, g_mnmpOpt->mask_level, i_chainsBatch.m_numReg[l_segOff+l_sg], i_chainsBatch.m_reg[l_segOff+l_sg], g_mnmpOpt->a * 2 + g_mnmpOpt->b, g_mnmpOpt->flag&MM_F_HARD_MLEVEL);
    mm_select_sub(NULL, g_mnmpOpt->pri_ratio, g_minimizer->k*2, g_mnmpOpt->best_n, &i_chainsBatch.m_numReg[l_segOff+l_sg], i_chainsBatch.m_reg[l_segOff+l_sg]);
    mm_set_sam_pri(i_chainsBatch.m_numReg[l_segOff+l_sg], i_chainsBatch.m_reg[l_segOff+l_sg]);
  }

  mm_set_mapq(NULL, i_chainsBatch.m_numReg[l_segOff+l_sg], i_chainsBatch.m_reg[l_segOff+l_sg], g_mnmpOpt->min_chain_score, g_mnmpOpt->a, i_repLen, i_isSrFlag);

  if (l_sg + 1 == l_numSeg && l_numSeg > 1) {
    mm_seg_free(NULL, l_numSeg, i_seg);
    if (n_segs == 2 && g_mnmpOpt->pe_ori >= 0 && (g_mnmpOpt->flag&MM_F_CIGAR)) {
      mm_pair(NULL, i_fragGap, g_mnmpOpt->pe_bonus, g_mnmpOpt->a * 2 + g_mnmpOpt->b, g_mnmpOpt->a, qlens, l_numReg, i_chainsBatch.m_reg[l_segOff+l_sg]);
    }
  }
}
