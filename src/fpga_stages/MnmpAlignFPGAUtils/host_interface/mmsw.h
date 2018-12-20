/*
 * =====================================================================================
 *
 *       Filename:  mmsw.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/26/2018 11:35:14 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiayi Sheng (jys), jysheng@falcon-computing.com
 *   Organization:  
 *
 * =====================================================================================
 */
#ifndef MMSW_H
#define MMSW_H
// #include "hls_stream.h"
// #include "ap_int.h"
// #include "ap_utils.h"
#include <stdint.h>

#define MAX_SEQ_LENGTH 256
#define MAX_ANKOR_NUM 2048
#define OUTPUT_UNIT_SIZE 1208
#define INPUT_UNIT_SIZE 33464
#define REGION_T_SIZE 80 // this is the mm_reg1_t size on the host side

#define KSW_NEG_INF -0x40000000
#define MM_PARENT_TMP_PRI (-2)
#define MM_SEED_LONG_JOIN  (1ULL<<40)
#define MAX_N_SEQ 128

#define KSW_EZ_SCORE_ONLY  0x01 // don't record alignment path/cigar
#define KSW_EZ_RIGHT       0x02 // right-align gaps
#define KSW_EZ_GENERIC_SC  0x04 // without this flag: match/mismatch only; last symbol is a wildcard
#define KSW_EZ_APPROX_MAX  0x08 // approximate max; this is faster with sse
#define KSW_EZ_APPROX_DROP 0x10 // approximate Z-drop; faster with sse
#define KSW_EZ_EXTZ_ONLY   0x40 // only perform extension
#define KSW_EZ_REV_CIGAR   0x80 // reverse CIGAR in the output
#define KSW_EZ_SPLICE_FOR  0x100
#define KSW_EZ_SPLICE_REV  0x200
#define KSW_EZ_SPLICE_FLANK 0x400



typedef struct {
	uint32_t capacity;                  // the capacity of cigar[]
	int32_t dp_score, dp_max, dp_max2;  // DP score; score of the max-scoring segment; score of the best alternate mappings
	//uint32_t n_ambi:30, trans_strand:2; // number of ambiguous bases; transcript strand: 0 for unknown, 1 for +, 2 for -
	uint32_t n_ambi;
    uint32_t trans_strand;
    uint32_t n_cigar;                   // number of cigar operations in cigar[]
    uint32_t cigar[MAX_SEQ_LENGTH];
} cigar_info;

typedef struct {
	uint32_t max:31, zdropped:1;
	int max_q, max_t;      // max extension coordinate
	int mqe, mqe_t;        // max score when reaching the end of query
	int mte, mte_q;        // max score when reaching the end of target
	int score;             // max score reaching both ends; may be KSW_NEG_INF
	int m_cigar, n_cigar;
	int reach_end;
    //uint32_t* cigar;
} extz_t;

typedef struct {
    uint64_t x;
    uint64_t y;
} anchor_t;

typedef struct {
	int32_t id;             // ID for internal uses (see also parent below)
	int32_t cnt;            // number of minimizers; if on the reverse strand
	int32_t rid;            // reference index; if this is an alignment from inversion rescue
	int32_t score;          // DP alignment score
	int32_t qs, qe, rs, re; // query start and end; reference start and end
	int32_t parent, subsc;  // parent==id if primary; best alternate mapping score
	int32_t as;             // offset in the a[] array (for internal uses only)
	int32_t mlen, blen;     // seeded exact match length; seeded alignment block length
	int32_t n_sub;          // number of suboptimal mappings
	int32_t score0;         // initial chaining score (before chain merging/spliting)
	//uint32_t mapq:8, split:2, rev:1, inv:1, sam_pri:1, proper_frag:1, pe_thru:1, seg_split:1, seg_id:8, split_inv:1, dummy:7;
	uint32_t mapq;
    uint32_t split;
    uint32_t rev;
    uint32_t inv;
    uint32_t sam_pri;
    uint32_t proper_frag;
    uint32_t pe_thru;
    uint32_t seg_split;
    uint32_t seg_id;
    uint32_t split_inv;
    uint32_t dummy;
	uint32_t hash;
	float div;
} region_t;

typedef union {
    float f;
    uint32_t d;
}float_union;

typedef struct {
	int seed;
	int sdust_thres; // score threshold for SDUST; 0 to disable
	int flag;        // see MM_F_* macros

	int bw;          // bandwidth
	int max_gap, max_gap_ref; // break a chain if there are no minimizers in a max_gap window
	int max_frag_len;
	int max_chain_skip;
	int min_cnt;         // min number of minimizers on each chain
	int min_chain_score; // min chaining score

	float mask_level;
	float pri_ratio;
	int best_n;      // top best_n chains are subjected to DP alignment

	int max_join_long, max_join_short;
	int min_join_flank_sc;
	float min_join_flank_ratio;

	int a, b, q, e, q2, e2; // matching score, mismatch, gap-open and gap-ext penalties
	int sc_ambi; // score when one or both bases are "N"
	int noncan;      // cost of non-canonical splicing sites
	int zdrop, zdrop_inv;   // break alignment if alignment score drops too fast along the diagonal
	int end_bonus;
	int min_dp_max;  // drop an alignment if the score of the max scoring segment is below this threshold
	int min_ksw_len;
	int anchor_ext_len, anchor_ext_shift;
	float max_clip_ratio; // drop an alignment if BOTH ends are clipped above this ratio

	int pe_ori, pe_bonus;

	float mid_occ_frac;  // only used by mm_mapopt_update(); see below
	int32_t min_mid_occ;
	int32_t mid_occ;     // ignore seeds with occurrences above this threshold
	int32_t max_occ;
	int mini_batch_size; // size of a batch of query bases to process in parallel
} opt_t;

typedef union {
    region_t orig;
    uint32_t serial[17]; 
} region_t_union;

typedef union {
    opt_t orig;
    uint32_t serial[41];
} opt_t_union;

typedef struct {
    int32_t data_size; // the entire size of current struct (taking n_a into consideration)
    int32_t qlen;
    int32_t n_a;
    uint32_t dp_score[2];
    region_t region[2];
    uint8_t qseq[2][MAX_SEQ_LENGTH];
    anchor_t a[MAX_ANKOR_NUM];
} align_input_fpga;

typedef struct {
    region_t_union region[2];
    cigar_info p;
} align_output_fpga;

typedef union {
    align_output_fpga out;
    uint32_t serial[OUTPUT_UNIT_SIZE / 4];
} align_output_union_fpga; 

extern "C" {
    void mmsw_kernel(const uint32_t* opt, int n_ref_seq, const uint64_t* ref_seqs, const uint32_t* ref, uint32_t* inputs, uint32_t* outputs, int batch_size);
}
#endif
