/*
 * =====================================================================================
 *
 *       Filename:  mmsw_baseline.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  11/13/2018 04:39:13 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiayi Sheng (jys), jysheng@falcon-computing.com
 *   Organization:  
 *
 * =====================================================================================
 */
#include "mmsw_baseline.h"
//#include "ksw2.h"
#include "kalloc.h"
#include <string.h>
#include <emmintrin.h>
#include <smmintrin.h>
using namespace std;

static inline void my_ksw_reset_extz(ksw_extz_t *ez)
{
  ez->max_q = ez->max_t = ez->mqe_t = ez->mte_q = -1;
  ez->max = 0, ez->score = ez->mqe = ez->mte = KSW_NEG_INF;
  ez->n_cigar = 0, ez->zdropped = 0, ez->reach_end = 0;
}

static inline int my_ksw_apply_zdrop(ksw_extz_t *ez, int is_rot, int32_t H, int a, int b, int zdrop, int8_t e, int idx)
{
  int r, t;
  if (is_rot) r = a, t = b;
  else r = a + b, t = a;
  //else r = a + b, t = a;
  if (H > (int32_t)ez->max) {
    ez->max = H, ez->max_t = t, ez->max_q = r - t;
  } else if (t >= ez->max_t && r - t >= ez->max_q) {
    int tl = t - ez->max_t, ql = (r - t) - ez->max_q, l;
    l = tl > ql? tl - ql : ql - tl;
    if (zdrop >= 0 && ez->max - H > zdrop + l * e) {
      ez->zdropped = 1;
      return 1;
    }
  }
  return 0;
}

static inline uint32_t *my_push_cigar(void *km, int *n_cigar, int *m_cigar, uint32_t *cigar, uint32_t op, int len)
{
  if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1]&0xf)) {
    if (*n_cigar > MAX_SEQ_LENGTH) {
      printf("n_cigar is large than 256, exiting\n");
      exit(-1);
    }
    if (*n_cigar == *m_cigar) {
      *m_cigar = *m_cigar? (*m_cigar)<<1 : 4;
      //cigar = (uint32_t*)krealloc(km, cigar, (*m_cigar) << 2);
    }
    cigar[(*n_cigar)++] = len<<4 | op;
  } else cigar[(*n_cigar)-1] += len<<4;
  return cigar;
}

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

typedef struct {
  int qlen, slen;
  uint8_t shift, mdiff, max, size;
  __m128i *qp, *H0, *H1, *E, *Hmax;
} kswq_t;

/**
 * Initialize the query data structure
 *
 * @param size   Number of bytes used to store a score; valid valures are 1 or 2
 * @param qlen   Length of the query sequence
 * @param query  Query sequence
 * @param m      Size of the alphabet
 * @param mat    Scoring matrix in a one-dimension array
 *
 * @return       Query data structure
 */
void *ksw_ll_qinit(void *km, int size, int qlen, const uint8_t *query, int m, const int8_t *mat)
{
  kswq_t *q;
  int slen, a, tmp, p;

  size = size > 1? 2 : 1;
  p = 8 * (3 - size); // # values per __m128i
  slen = (qlen + p - 1) / p; // segmented length
  q = (kswq_t*)kmalloc(km, sizeof(kswq_t) + 256 + 16 * slen * (m + 4)); // a single block of memory
  q->qp = (__m128i*)(((size_t)q + sizeof(kswq_t) + 15) >> 4 << 4); // align memory
  q->H0 = q->qp + slen * m;
  q->H1 = q->H0 + slen;
  q->E  = q->H1 + slen;
  q->Hmax = q->E + slen;
  q->slen = slen; q->qlen = qlen; q->size = size;
  // compute shift
  tmp = m * m;
  for (a = 0, q->shift = 127, q->mdiff = 0; a < tmp; ++a) { // find the minimum and maximum score
    if (mat[a] < (int8_t)q->shift) q->shift = mat[a];
    if (mat[a] > (int8_t)q->mdiff) q->mdiff = mat[a];
  }
  q->max = q->mdiff;
  q->shift = 256 - q->shift; // NB: q->shift is uint8_t
  q->mdiff += q->shift; // this is the difference between the min and max scores
  // An example: p=8, qlen=19, slen=3 and segmentation:
  //  {{0,3,6,9,12,15,18,-1},{1,4,7,10,13,16,-1,-1},{2,5,8,11,14,17,-1,-1}}
  if (size == 1) {
    int8_t *t = (int8_t*)q->qp;
    for (a = 0; a < m; ++a) {
      int i, k, nlen = slen * p;
      const int8_t *ma = mat + a * m;
      for (i = 0; i < slen; ++i)
        for (k = i; k < nlen; k += slen) // p iterations
          *t++ = (k >= qlen? 0 : ma[query[k]]) + q->shift;
    }
  } else {
    int16_t *t = (int16_t*)q->qp;
    for (a = 0; a < m; ++a) {
      int i, k, nlen = slen * p;
      const int8_t *ma = mat + a * m;
      for (i = 0; i < slen; ++i)
        for (k = i; k < nlen; k += slen) // p iterations
          *t++ = (k >= qlen? 0 : ma[query[k]]);
    }
  }
  return q;
}

int ksw_ll_i16(void *q_, int tlen, const uint8_t *target, int _gapo, int _gape, int *qe, int *te)
{
  kswq_t *q = (kswq_t*)q_;
  int slen, i, gmax = 0, qlen8;
  __m128i zero, gapoe, gape, *H0, *H1, *E, *Hmax;
  uint16_t *H8;

#define __max_8(ret, xx) do { \
    (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 8)); \
    (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 4)); \
    (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 2)); \
    (ret) = _mm_extract_epi16((xx), 0); \
  } while (0)

  // initialization
  *qe = *te = -1;
  zero = _mm_set1_epi32(0);
  gapoe = _mm_set1_epi16(_gapo + _gape);
  gape = _mm_set1_epi16(_gape);
  H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;
  slen = q->slen, qlen8 = slen * 8;
  memset(E,    0, slen * sizeof(__m128i));
  memset(H0,   0, slen * sizeof(__m128i));
  memset(Hmax, 0, slen * sizeof(__m128i));
  // the core loop
  for (i = 0; i < tlen; ++i) {
    int j, k, imax;
    __m128i e, h, f = zero, max = zero, *S = q->qp + target[i] * slen; // s is the 1st score vector
    h = _mm_load_si128(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
    h = _mm_slli_si128(h, 2);
    for (j = 0; LIKELY(j < slen); ++j) {
      h = _mm_adds_epi16(h, *S++);
      e = _mm_load_si128(E + j);
      h = _mm_max_epi16(h, e);
      h = _mm_max_epi16(h, f);
      max = _mm_max_epi16(max, h);
      _mm_store_si128(H1 + j, h);
      h = _mm_subs_epu16(h, gapoe);
      e = _mm_subs_epu16(e, gape);
      e = _mm_max_epi16(e, h);
      _mm_store_si128(E + j, e);
      f = _mm_subs_epu16(f, gape);
      f = _mm_max_epi16(f, h);
      h = _mm_load_si128(H0 + j);
    }
    for (k = 0; LIKELY(k < 8); ++k) {
      f = _mm_slli_si128(f, 2);
      for (j = 0; LIKELY(j < slen); ++j) {
        h = _mm_load_si128(H1 + j);
        h = _mm_max_epi16(h, f);
        _mm_store_si128(H1 + j, h);
        h = _mm_subs_epu16(h, gapoe);
        f = _mm_subs_epu16(f, gape);
        if(UNLIKELY(!_mm_movemask_epi8(_mm_cmpgt_epi16(f, h)))) goto end_loop_i16;
      }
    }
end_loop_i16:
    __max_8(imax, max);
    if (imax >= gmax) {
      gmax = imax; *te = i;
      memcpy(Hmax, H1, slen * sizeof(__m128i));
    }
    S = H1; H1 = H0; H0 = S;
  }
  for (i = 0, H8 = (uint16_t*)Hmax; i < qlen8; ++i)
    if ((int)H8[i] == gmax) *qe = i / 8 + i % 8 * slen;
  return gmax;
}

int ksw_ll(void *q_, int tlen, const uint8_t *target, int _gapo, int _gape, int *qe, int *te)
{
  kswq_t *q = (kswq_t*)q_;
  int slen, i, gmax = 0, qlen8;
  __m128i zero, gapoe, gape, *H0, *H1, *E, *Hmax;
  uint16_t *H8;

#define __max_8(ret, xx) do { \
    (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 8)); \
    (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 4)); \
    (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 2)); \
    (ret) = _mm_extract_epi16((xx), 0); \
  } while (0)

  // initialization
  *qe = *te = -1;
  zero = _mm_set1_epi32(0);
  gapoe = _mm_set1_epi16(_gapo + _gape);
  gape = _mm_set1_epi16(_gape);
  H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;
  slen = q->slen, qlen8 = slen * 8;
  memset(E,    0, slen * sizeof(__m128i));
  memset(H0,   0, slen * sizeof(__m128i));
  memset(Hmax, 0, slen * sizeof(__m128i));
  // the core loop
  for (i = 0; i < tlen; ++i) {
    int j, k, imax;
    __m128i e, h, f = zero, max = zero, *S = q->qp + target[i] * slen; // s is the 1st score vector
    h = _mm_load_si128(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
    h = _mm_slli_si128(h, 2);
    for (j = 0; LIKELY(j < slen); ++j) {
      h = _mm_adds_epi16(h, *S++);
      e = _mm_load_si128(E + j);
      h = _mm_max_epi16(h, e);
      h = _mm_max_epi16(h, f);
      max = _mm_max_epi16(max, h);
      _mm_store_si128(H1 + j, h);
      h = _mm_subs_epu16(h, gapoe);
      e = _mm_subs_epu16(e, gape);
      e = _mm_max_epi16(e, h);
      _mm_store_si128(E + j, e);
      f = _mm_subs_epu16(f, gape);
      f = _mm_max_epi16(f, h);
      h = _mm_load_si128(H0 + j);
    }
    for (k = 0; LIKELY(k < 8); ++k) {
      f = _mm_slli_si128(f, 2);
      for (j = 0; LIKELY(j < slen); ++j) {
        h = _mm_load_si128(H1 + j);
        h = _mm_max_epi16(h, f);
        _mm_store_si128(H1 + j, h);
        h = _mm_subs_epu16(h, gapoe);
        f = _mm_subs_epu16(f, gape);
        if(UNLIKELY(!_mm_movemask_epi8(_mm_cmpgt_epi16(f, h)))) goto end_loop_i16;
      }
    }
end_loop_i16:
    __max_8(imax, max);
    if (imax >= gmax) {
      gmax = imax; *te = i;
      memcpy(Hmax, H1, slen * sizeof(__m128i));
    }
    S = H1; H1 = H0; H0 = S;
  }
  for (i = 0, H8 = (uint16_t*)Hmax; i < qlen8; ++i)
    if ((int)H8[i] == gmax) *qe = i / 8 + i % 8 * slen;
  return gmax;
}

#define mm_seq4_get(s, i)    ((s)[(i)>>3] >> (((i)&7)<<2) & 0xf)

int my_mm_idx_getseq(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq)
{
  uint64_t i, st1, en1;
  if (rid >= mi->n_seq || st >= mi->seq[rid].len) return -1;
  if (en > mi->seq[rid].len) en = mi->seq[rid].len;
  st1 = mi->seq[rid].offset + st;
  en1 = mi->seq[rid].offset + en;
  for (i = st1; i < en1; ++i)
    seq[i - st1] = mm_seq4_get(mi->S, i);
  return en - st;
}

static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b, int8_t sc_ambi)
{
  //generating matrix that can be used for look up with match score and mismatch penalty
  //a is the match score, b is the mismatch penalty
  //ACGT and N corresponding to 0123 and 5
  //sc_ambi is the score of ambiguise. the fifth colume and fifth row are populated with sc_ambi
  int i, j;
  a = a < 0? -a : a;
  b = b > 0? -b : b;
  sc_ambi = sc_ambi > 0? -sc_ambi : sc_ambi;
  for (i = 0; i < m - 1; ++i) {
    for (j = 0; j < m - 1; ++j)
      mat[i * m + j] = i == j? a : b;
    mat[i * m + m - 1] = sc_ambi;
  }
  for (j = 0; j < m; ++j)
    mat[(m - 1) * m + j] = sc_ambi;
}

static inline void mm_seq_rev(uint32_t len, uint8_t *seq)
{
  uint32_t i;
  uint8_t t;
  for (i = 0; i < len>>1; ++i)
    t = seq[i], seq[i] = seq[len - 1 - i], seq[len - 1 - i] = t;
}

static void mm_fix_cigar(mm_reg1_t *r, const uint8_t *qseq, const uint8_t *tseq, int *qshift, int *tshift)
{
  mm_extra_t *p = r->p;
  int32_t toff = 0, qoff = 0, to_shrink = 0;
  uint32_t k;
  *qshift = *tshift = 0;
  if (p->n_cigar <= 1) return;
  for (k = 0; k < p->n_cigar; ++k) { // indel left alignment
    uint32_t op = p->cigar[k]&0xf, len = p->cigar[k]>>4;
    if (len == 0) to_shrink = 1;
    if (op == 0) {
      toff += len, qoff += len;
    } else if (op == 1 || op == 2) { // insertion or deletion
      if (k > 0 && k < p->n_cigar - 1 && (p->cigar[k-1]&0xf) == 0 && (p->cigar[k+1]&0xf) == 0) {
        int l, prev_len = p->cigar[k-1] >> 4;
        if (op == 1) {
          for (l = 0; l < prev_len; ++l)
            if (qseq[qoff - 1 - l] != qseq[qoff + len - 1 - l])
              break;
        } else {
          for (l = 0; l < prev_len; ++l)
            if (tseq[toff - 1 - l] != tseq[toff + len - 1 - l])
              break;
        }
        if (l > 0)
          p->cigar[k-1] -= l<<4, p->cigar[k+1] += l<<4, qoff -= l, toff -= l;
        if (l == prev_len) to_shrink = 1;
      }
      if (op == 1) qoff += len;
      else toff += len;
    } else if (op == 3) {
      toff += len;
    }
  }
  //assert(qoff == r->qe - r->qs && toff == r->re - r->rs);
  if (to_shrink) { // squeeze out zero-length operations
    int32_t l = 0;
    for (k = 0; k < p->n_cigar; ++k) // squeeze out zero-length operations
      if (p->cigar[k]>>4 != 0)
        p->cigar[l++] = p->cigar[k];
    p->n_cigar = l;
    for (k = l = 0; k < p->n_cigar; ++k) // merge two adjacent operations if they are the same
      if (k == p->n_cigar - 1 || (p->cigar[k]&0xf) != (p->cigar[k+1]&0xf))
        p->cigar[l++] = p->cigar[k];
      else p->cigar[k+1] += p->cigar[k]>>4<<4; // add length to the next CIGAR operator
    p->n_cigar = l;
  }
  if ((p->cigar[0]&0xf) == 1 || (p->cigar[0]&0xf) == 2) { // get rid of leading I or D
    int32_t l = p->cigar[0] >> 4;
    if ((p->cigar[0]&0xf) == 1) {
      if (r->rev) r->qe -= l;
      else r->qs += l;
      *qshift = l;
    } else r->rs += l, *tshift = l;
    --p->n_cigar;
    memmove(p->cigar, p->cigar + 1, p->n_cigar * 4);
  }
}

static inline void mm_cal_fuzzy_len(mm_reg1_t *r, const mm128_t *a)
{
  int i;
  r->mlen = r->blen = 0;
  if (r->cnt <= 0) return;
  r->mlen = r->blen = a[r->as].y>>32&0xff;
  for (i = r->as + 1; i < r->as + r->cnt; ++i) {
    int span = a[i].y>>32&0xff;
    int tl = (int32_t)a[i].x - (int32_t)a[i-1].x;
    int ql = (int32_t)a[i].y - (int32_t)a[i-1].y;
    r->blen += tl > ql? tl : ql;
    r->mlen += tl > span && ql > span? span : tl < ql? tl : ql;
  }
}

static inline void mm_reg_set_coor(mm_reg1_t *r, int32_t qlen, const mm128_t *a)
{ // NB: r->as and r->cnt MUST BE set correctly for this function to work
  int32_t k = r->as, q_span = (int32_t)(a[k].y>>32&0xff);
  r->rev = a[k].x>>63;
  r->rid = a[k].x<<1>>33;
  r->rs = (int32_t)a[k].x + 1 > q_span? (int32_t)a[k].x + 1 - q_span : 0; // NB: target span may be shorter, so this test is necessary
  r->re = (int32_t)a[k + r->cnt - 1].x + 1;
  if (!r->rev) {
    r->qs = (int32_t)a[k].y + 1 - q_span;
    r->qe = (int32_t)a[k + r->cnt - 1].y + 1;
  } else {
    r->qs = qlen - ((int32_t)a[k + r->cnt - 1].y + 1);
    r->qe = qlen - ((int32_t)a[k].y + 1 - q_span);
  }
  mm_cal_fuzzy_len(r, a);
}

void my_mm_split_reg(mm_reg1_t *r, mm_reg1_t *r2, int n, int qlen, mm128_t *a)
{
  if (n <= 0 || n >= r->cnt) return;
  *r2 = *r;
  r2->id = -1;
  r2->sam_pri = 0;
  r2->p = 0;
  r2->split_inv = 0;
  r2->cnt = r->cnt - n;
  r2->score = (int32_t)(r->score * ((float)r2->cnt / r->cnt) + .499);
  r2->as = r->as + n;
  if (r->parent == r->id) r2->parent = MM_PARENT_TMP_PRI;
  mm_reg_set_coor(r2, qlen, a);
  r->cnt -= r2->cnt;
  r->score -= r2->score;
  mm_reg_set_coor(r, qlen, a);
  r->split |= 1, r2->split |= 2;
}

static void mm_max_stretch(const mm_reg1_t *r, const mm128_t *a, int32_t *as, int32_t *cnt)
{
  int32_t i, score, max_score, len, max_i, max_len;

  *as = r->as, *cnt = r->cnt;
  if (r->cnt < 2) return;

  max_score = -1, max_i = -1, max_len = 0;
  score = a[r->as].y >> 32 & 0xff, len = 1;
  for (i = r->as + 1; i < r->as + r->cnt; ++i) {
    int32_t lq, lr, q_span;
    q_span = a[i].y >> 32 & 0xff;
    lr = (int32_t)a[i].x - (int32_t)a[i-1].x;
    lq = (int32_t)a[i].y - (int32_t)a[i-1].y;
    if (lq == lr) {
      score += lq < q_span? lq : q_span;
      ++len;
    } else {
      if (score > max_score)
        max_score = score, max_len = len, max_i = i - len;
      score = q_span, len = 1;
    }
  }
  if (score > max_score)
    max_score = score, max_len = len, max_i = i - len;
  *as = max_i, *cnt = max_len;
}

// In the backtrack matrix, value p[] has the following structure:
//   bit 0-2: which type gets the max - 0 for H, 1 for E, 2 for F, 3 for \tilde{E} and 4 for \tilde{F}
//   bit 3/0x08: 1 if a continuation on the E state (bit 5/0x20 for a continuation on \tilde{E})
//   bit 4/0x10: 1 if a continuation on the F state (bit 6/0x40 for a continuation on \tilde{F})
static inline void my_ksw_backtrack(void *km, int is_rot, int is_rev, int min_intron_len, const uint8_t *p, const int *off, const int *off_end, int n_col, int i0, int j0,
                 int *m_cigar_, int *n_cigar_, uint32_t **cigar_)
{ // p[] - lower 3 bits: which type gets the max; bit
  int n_cigar = 0, m_cigar = *m_cigar_, i = i0, j = j0, r, state = 0;
  uint32_t *cigar = *cigar_, tmp;
  while (i >= 0 && j >= 0) { // at the beginning of the loop, _state_ tells us which state to check
    int force_state = -1;
    if (is_rot) {
      r = i + j;
      if (i < off[r]) force_state = 2;
      if (off_end && i > off_end[r]) force_state = 1;
      tmp = force_state < 0? p[(size_t)r * n_col + i - off[r]] : 0;
    } else {
      if (j < off[i]) force_state = 2;
      if (off_end && j > off_end[i]) force_state = 1;
      tmp = force_state < 0? p[(size_t)i * n_col + j - off[i]] : 0;
    }
    if (state == 0) state = tmp & 7; // if requesting the H state, find state one maximizes it.
    else if (!(tmp >> (state + 2) & 1)) state = 0; // if requesting other states, _state_ stays the same if it is a continuation; otherwise, set to H
    if (state == 0) state = tmp & 7; // TODO: probably this line can be merged into the "else if" line right above; not 100% sure
    if (force_state >= 0) state = force_state;
    if (state == 0) cigar = my_push_cigar(km, &n_cigar, &m_cigar, cigar, 0, 1), --i, --j; // match
    else if (state == 1 || (state == 3 && min_intron_len <= 0)) cigar = my_push_cigar(km, &n_cigar, &m_cigar, cigar, 2, 1), --i; // deletion
    else if (state == 3 && min_intron_len > 0) cigar = my_push_cigar(km, &n_cigar, &m_cigar, cigar, 3, 1), --i; // intron
    else cigar = my_push_cigar(km, &n_cigar, &m_cigar, cigar, 1, 1), --j; // insertion
  }
  if (i >= 0) cigar = my_push_cigar(km, &n_cigar, &m_cigar, cigar, min_intron_len > 0 && i >= min_intron_len? 3 : 2, i + 1); // first deletion
  if (j >= 0) cigar = my_push_cigar(km, &n_cigar, &m_cigar, cigar, 1, j + 1); // first insertion
  if (!is_rev)
    for (i = 0; i < n_cigar>>1; ++i) // reverse CIGAR
      tmp = cigar[i], cigar[i] = cigar[n_cigar-1-i], cigar[n_cigar-1-i] = tmp;
  *m_cigar_ = m_cigar, *n_cigar_ = n_cigar, *cigar_ = cigar;
}

typedef struct { int32_t h, e, e2; } eh_t;

void ksw_extd(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
        int8_t gapo, int8_t gape, int8_t gapo2, int8_t gape2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez, int idx)
{
  eh_t *eh;
  int8_t *qp; // query profile
  int32_t i, j, k, max_j = 0, gapoe = gapo + gape, gapoe2 = gapo2 + gape2, n_col, *off = 0, with_cigar = !(flag&KSW_EZ_SCORE_ONLY);
  uint8_t *z = 0; // backtrack matrix; in each cell: f<<4|e<<2|h; in principle, we can halve the memory, but backtrack will be more complex

  my_ksw_reset_extz(ez);

  // allocate memory
  if (w < 0) w = tlen > qlen? tlen : qlen;
  n_col = qlen < 2*w+1? qlen : 2*w+1; // maximum #columns of the backtrack matrix
  qp = (int8_t*)kmalloc(km, qlen * m);
  eh = (eh_t*)kcalloc(km, qlen + 1, sizeof(eh_t));
  if (with_cigar) {
    z = (uint8_t*)kmalloc(km, (size_t)n_col * tlen);
    off = (int32_t*)kcalloc(km, tlen, 4);
  }

  // generate the query profile
  for (k = i = 0; k < m; ++k) {
    const int8_t *p = &mat[k * m];
    for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
  }

  // fill the first row
  eh[0].h = 0, eh[0].e = -gapoe - gapoe, eh[0].e2 = -gapoe2 - gapoe2;
  for (j = 1; j <= qlen && j <= w; ++j) {
    int tmp;
    eh[j].h = -(gapo + gape * j) > -(gapo2 + gape2 * j)? -(gapo + gape * j) : -(gapo2 + gape2 * j);
    tmp = -(gapoe + gape * j) > -(gapoe2 + gape2 * j)? -(gapoe + gape * j) : -(gapoe2 + gape2 * j);
    eh[j].e = tmp - gapoe;
    eh[j].e2 = tmp - gapoe2;
  }
  for (; j <= qlen; ++j) eh[j].h = eh[j].e = eh[j].e2 = KSW_NEG_INF; // everything is -inf outside the band

  for (int r = 0; r < qlen + tlen - 1; ++r) {
    int st = 0, en = tlen - 1;
    if (st < r - qlen + 1) st = r - qlen + 1;
    if (en > r) en = r;
    if (st < (r-w+1)>>1) st = (r-w+1)>>1; // take the ceil
    if (en > (r+w)>>1) en = (r+w)>>1; // take the floor
    if (st > en) {
      ez->zdropped = 1;
      break;
    }
  }
  // DP loop
  if (!ez->zdropped) {
  for (i = 0; i < tlen; ++i) { // target sequence is in the outer loop
    
    int32_t f, f2, h1, st, en, max = KSW_NEG_INF, tmp;
    int8_t *q = &qp[target[i] * qlen];
    st = i > w? i - w : 0;
    en = i + w < qlen - 1? i + w : qlen - 1;
    tmp = -(gapoe + gape * i) > -(gapoe2 + gape2 * i)? -(gapoe + gape * i) : -(gapoe2 + gape2 * i);
    h1 = st > 0? KSW_NEG_INF : tmp;
    f  = st > 0? KSW_NEG_INF : tmp - gapoe;
    f2 = st > 0? KSW_NEG_INF : tmp - gapoe2;
    if (!with_cigar) {
      for (j = st; j <= en; ++j) {
   
        eh_t *p = &eh[j];
        int32_t h = p->h, h2, e = p->e, e2 = p->e2;
        p->h = h1;
        h += q[j];
        h = h >= e?  h : e;
        h = h >= f?  h : f;
        h = h >= e2? h : e2;
        h = h >= f2? h : f2;
        h1 = h;
        max_j = max > h? max_j : j;
        max   = max > h? max   : h;
        h -= gapoe;
        e -= gape;
        e  = e > h? e : h;
        p->e = e;
        f -= gape;
        f  = f > h? f : h;
        h2 = h1 - gapoe2;
        e2-= gape2;
        e2 = e2 > h2? e2 : h2;
        p->e2 = e2;
        f2-= gape2;
        f2 = f2 > h2? f2 : h2;
      }
    } else if (!(flag&KSW_EZ_RIGHT)) {
      uint8_t *zi = &z[(long)i * n_col];
      off[i] = st;
      for (j = st; j <= en; ++j) {

        eh_t *p = &eh[j];
        int32_t h = p->h, h2, e = p->e, e2 = p->e2;
        uint8_t d; // direction
        p->h = h1;
        h += q[j];
        d = h >= e?  0 : 1;
        h = h >= e?  h : e;
        d = h >= f?  d : 2;
        h = h >= f?  h : f;
        d = h >= e2? d : 3;
        h = h >= e2? h : e2;
        d = h >= f2? d : 4;
        h = h >= f2? h : f2;
        h1 = h;
        max_j = max > h? max_j : j;
        max   = max > h? max   : h;
        h -= gapoe;
        e -= gape;
        d |= e > h? 1<<3 : 0;
        e  = e > h? e    : h;
        p->e = e;
        f -= gape;
        d |= f > h? 1<<4 : 0; // if we want to halve the memory, use one bit only, instead of two
        f  = f > h? f    : h;
        h2 = h1 - gapoe2;
        e2-= gape2;
        d |= e2 > h2? 1<<5 : 0;
        e2 = e2 > h2? e2 : h2;
        p->e2 = e2;
        f2-= gape2;
        d |= f2 > h2? 1<<6 : 0;
        f2 = f2 > h2? f2 : h2;
        zi[j - st] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
      }
    } else {
      uint8_t *zi = &z[(long)i * n_col];
      off[i] = st;
      for (j = st; j <= en; ++j) {

        eh_t *p = &eh[j];
        int32_t h = p->h, h2, e = p->e, e2 = p->e2;
        uint8_t d; // direction
        p->h = h1;
        h += q[j];
        d = h > e?  0 : 1;
        h = h > e?  h : e;
        d = h > f?  d : 2;
        h = h > f?  h : f;
        d = h > e2? d : 3;
        h = h > e2? h : e2;
        d = h > f2? d : 4;
        h = h > f2? h : f2;
        h1 = h;
        max_j = max > h? max_j : j;
        max   = max > h? max   : h;
        h -= gapoe;
        e -= gape;
        d |= e >= h? 1<<3 : 0;
        e  = e >= h? e    : h;
        p->e = e;
        f -= gape;
        d |= f >= h? 1<<4 : 0; // if we want to halve the memory, use one bit only, instead of two
        f  = f >= h? f    : h;
        h2 = h1 - gapoe2;
        e2-= gape2;
        d |= e2 >= h2? 1<<5 : 0;
        e2 = e2 >= h2? e2 : h2;
        p->e2 = e2;
        f2-= gape2;
        d |= f2 >= h2? 1<<6 : 0;
        f2 = f2 >= h2? f2 : h2;
        zi[j - st] = d; // z[i,j] keeps h for the current cell and e/f for the next cell
      }
    }
    eh[j].h = h1, eh[j].e = KSW_NEG_INF;
    // update ez
    if (en == qlen - 1 && eh[qlen].h > ez->mqe)
      ez->mqe = eh[qlen].h, ez->mqe_t = i;
    if (i == tlen - 1)
      ez->mte = max, ez->mte_q = max_j;
    if (my_ksw_apply_zdrop(ez, 0, max, i, max_j, zdrop, gape2, idx)) break;
    if (i == tlen - 1 && en == qlen - 1)
      ez->score = eh[qlen].h;
  }
  }
  kfree(km, qp); kfree(km, eh);
  if (with_cigar) {
    int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
    if (!ez->zdropped && !(flag&KSW_EZ_EXTZ_ONLY))
      my_ksw_backtrack(km, 0, rev_cigar, 0, z, off, 0, n_col, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
    else if (!ez->zdropped && (flag&KSW_EZ_EXTZ_ONLY) && ez->mqe + end_bonus > (int)ez->max) {
      ez->reach_end = 1;
      my_ksw_backtrack(km, 0, rev_cigar, 0, z, off, 0, n_col, ez->mqe_t, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
    }
    else if (ez->max_t >= 0 && ez->max_q >= 0)
      my_ksw_backtrack(km, 0, rev_cigar, 0, z, off, 0, n_col, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
    kfree(km, z); kfree(km, off);
  }

}

void ksw_gg2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez, int idx)
{
  int r, t, qe = q + e, qe2 = q2 + e2, n_col_, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en, wl, wr, max_sc, min_sc, long_thres, long_diff;
  int with_cigar = !(flag&KSW_EZ_SCORE_ONLY), approx_max = !!(flag&KSW_EZ_APPROX_MAX);
  int8_t *u, *v, *x, *y, *x2, *y2, *s;
  int32_t *H = 0, H0 = 0, last_H0_t = 0;
  uint8_t *qr, *sf, *mem, *mem2 = 0;
  uint8_t *p = 0;

  my_ksw_reset_extz(ez);
  if (m <= 1 || qlen <= 0 || tlen <= 0) return;
  
  if (q2 + e2 < q + e) t = q, q = q2, q2 = t, t = e, e = e2, e2 = t; // make sure q+e no larger than q2+e2
  
  if (w < 0) w = tlen > qlen? tlen : qlen;
  wl = wr = w;
  int n_col = qlen < tlen? qlen : tlen;
  n_col = (n_col < w + 1)? n_col : w + 1;
  for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
    max_sc = max_sc > mat[t]? max_sc : mat[t];
    min_sc = min_sc < mat[t]? min_sc : mat[t];
  }
  if (-min_sc > 2 * (q + e)) return; // otherwise, we won't see any mismatches

  long_thres = e != e2? (q2 - q) / (e - e2) - 1 : 0;

  if (q2 + e2 + long_thres * e2 > q + e + long_thres * e)
    ++long_thres;
  long_diff = long_thres * (e - e2) - (q2 - q) - e2;

  mem = (uint8_t*)kcalloc(km, (tlen + 1) * 8 + qlen + 2, 1);
  u = (int8_t*)mem;
  v = u + tlen + 1, x = v + tlen + 1, y = x + tlen + 1, x2 = y + tlen + 1, y2 = x2 + tlen + 1;
  s = y2 + tlen + 1, sf = (uint8_t*)(s + tlen + 1), qr = sf + tlen + 1;
  memset(u,  -q  - e,  tlen + 1);
  memset(v,  -q  - e,  tlen + 1);
  memset(x,  -q  - e,  tlen + 1);
  memset(y,  -q  - e,  tlen + 1);
  memset(x2, -q2 - e2, tlen + 1);
  memset(y2, -q2 - e2, tlen + 1);

  if (!approx_max) {
    H = (int32_t*)kmalloc(km, (tlen + 1) * 4);
    for (t = 0; t < tlen + 1; ++t) H[t] = KSW_NEG_INF;
  }

  if (with_cigar) {
    p = (uint8_t*)kcalloc(km, (qlen + tlen) * n_col, 1);
    off = (int*)kmalloc(km, (qlen + tlen - 1) * sizeof(int) * 2);
    off_end = off + qlen + tlen - 1;
  }

  for (t = 0; t < qlen; ++t)
    qr[t] = query[qlen - 1 - t];
  
  memcpy(sf, target, tlen);

  for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r) {
    int st = 0, en = tlen - 1, st0, en0, st_, en_;
    int8_t x1, x21, v1;
    uint8_t *qrr = qr + (qlen - 1 - r);
    // find the boundaries
    if (st < r - qlen + 1) st = r - qlen + 1;
    if (en > r) en = r;
    if (st < (r-w+1)>>1) st = (r-w+1)>>1; // take the ceil
    if (en > (r+w)>>1) en = (r+w)>>1; // take the floor
    if (st > en) {
      printf("zdropped in st > en\n");
      ez->zdropped = 1;
      break;
    }
    st0 = st, en0 = en;
    
    // set boundary conditions
    if (st > 0) {
      if (st - 1 >= last_st && st - 1 <= last_en) {
        x1 = x[st - 1], x21 = x2[st - 1], v1 = v[st - 1]; // (r-1,s-1) calculated in the last round
      } else {
        x1 = -q - e, x21 = -q2 - e2;
        v1 = -q - e;
      }
    } else {
      x1 = -q - e, x21 = -q2 - e2;
      v1 = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
    }
    if (en >= r) {
      ((int8_t*)y)[r] = -q - e, ((int8_t*)y2)[r] = -q2 - e2;
      u[r] = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
    }
    
    
    // loop fission: set scores first
    if (!(flag & KSW_EZ_GENERIC_SC)) {
      uint8_t sc_mch, sc_mis, sc_N;
      sc_mch = mat[0]; sc_mis = mat[1]; sc_N = mat[m * m - 1];
      for (t = st; t <= en; ++t) {
        uint8_t tmp = (sf[t] == qrr[t]) ? sc_mch : sc_mis;
        if (sf[t] == m - 1 || qrr[t] == m - 1) {
          tmp = sc_N;
        }
        s[t] = tmp;
      }
    } else {
      for (t = st; t <= en; ++t)
        s[t] = mat[sf[t] * m + qrr[t]];
    }
    
    // core loop
    if (!with_cigar) {
      for (t = st; t <= en; ++t) {
        //dp_code_block1_start
        int8_t z = s[t];
        int8_t a = x1   + v1;
        int8_t b = y[t] + u[t];
        int8_t a2 = x21 + v1;
        int8_t b2 = y2[t] + u[t];
        //dp_code_block1_end
        z = a > z? a : z;
        z = b > z? b : z;
        z = a2 > z? a2 : z;
        z = b2 > z? b2 : z;
        z = mat[0] > z? z : mat[0];
        //dp_code_block2_start
        int8_t u1;
        u1 = u[t];
        u[t] = z - v1;
        v1 = v[t];
        v[t] = z - u1;
        int8_t tmp = z - q;
        a -= tmp;
        b -= tmp;
        tmp = z - q2;
        a2 -= tmp;
        b2 -= tmp;
        //dp_code_block2_end
        x1 = x[t];
        x21 = x2[t];
        x[t] = (a > 0 ? a : 0) - qe;
        y[t] = (b > 0 ? b : 0) - qe;
        x2[t] = (a2 > 0 ? a2 : 0) - qe2;
        y2[t] = (b2 > 0 ? b2 : 0) - qe2;
      }
    } else if (!(flag&KSW_EZ_RIGHT)) { // gap left-alignment
      uint8_t *pr = p + r * n_col - st;
      off[r] = st, off_end[r] = en;
      for (t = st; t <= en; ++t) {
        //dp_code_block1_start
        int8_t z = s[t];
        int8_t a = x1   + v1;
        int8_t b = y[t] + u[t];
        int8_t a2 = x21 + v1;
        int8_t b2 = y2[t] + u[t];
        //dp_code_block1_end
        uint8_t d;
        d = a > z? 1 : 0;
        z = a > z? a : z;
        d = b > z? 2 : d;
        z = b > z? b : z;
        d = a2 > z ? 3 : d;
        z = a2 > z ? a2 : z;
        d = b2 > z ? 4 : d;
        z = b2 > z ? b2 : z;
        z = z < mat[0] ? z : mat[0];
        //dp_code_block2_start
        int8_t u1;
        u1 = u[t];
        u[t] = z - v1;
        v1 = v[t];
        v[t] = z - u1;
        int8_t tmp = z - q;
        a -= tmp;
        b -= tmp;
        tmp = z - q2;
        a2 -= tmp;
        b2 -= tmp;
        //dp_code_block2_end
        x1 = x[t];
        x21 = x2[t];
        d   |= a > 0? 0x08 : 0;
        x[t] = (a > 0 ? a : 0) - qe;
        d   |= b > 0? 0x10 : 0;
        y[t] = (b > 0 ? b : 0) - qe;
        d   |= a2 > 0? 0x20 : 0;
        x2[t] = (a2 > 0 ? a2 : 0) - qe2;
        d   |= b2 > 0? 0x40 : 0;
        y2[t] = (b2 > 0 ? b2 : 0) - qe2;
        pr[t] = d;
      }
    } else {
      uint8_t *pr = p + r * n_col - st;
      off[r] = st, off_end[r] = en;
      for (t = st; t <= en; ++t) {
        //dp_code_block1_start
        int8_t z = s[t];
        int8_t a = x1   + v1;
        int8_t b = y[t] + u[t];
        int8_t a2 = x21 + v1;
        int8_t b2 = y2[t] + u[t];
        //dp_code_block1_end
        uint8_t d;
        d = z > a? 0 : 1;
        z = a > z? a : z;
        d = z > b? d : 2;
        z = b > z? b : z;
        d = z > a2 ? d : 3;
        z = a2 > z ? a2 : z;
        d = z > b2 ? d : 4;
        z = b2 > z ? b2 : z;
        z = z < mat[0] ? z : mat[0];
        //dp_code_block2_start
        int8_t u1;
        u1 = u[t];
        u[t] = z - v1;
        v1 = v[t];
        v[t] = z - u1;
        int8_t tmp = z - q;
        a -= tmp;
        b -= tmp;
        tmp = z - q2;
        a2 -= tmp;
        b2 -= tmp;
        //dp_code_block2_end
        x1 = x[t];
        x21 = x2[t];
        d   |= a >= 0? 0x08 : 0;
        x[t] = (a >= 0 ? a : 0) - qe;
        d   |= b >= 0? 0x10 : 0;
        y[t] = (b >= 0 ? b : 0) - qe;
        d   |= a2 >= 0? 0x20 : 0;
        x2[t] = (a2 >= 0 ? a2 : 0) - qe2;
        d   |= b2 >= 0? 0x40 : 0;
        y2[t] = (b2 >= 0 ? b2 : 0) - qe2;
        pr[t] = d;
      }
    }
    if (!approx_max) { // find the exact max with a 32-bit score array
      int32_t max_H, max_t;
      // compute H[], max_H and max_t
      if (r > 0) {
        max_H = H[en] = en > 0? H[en-1] + u[en] : H[en] + v[en]; // special casing the last element
        max_t = en;
        for (t = st; t < en; t++) {
          H[t] += v[t];
          if (H[t] > max_H) {
            max_H = H[t];
            max_t = t;
          }
        }
      } else {
        H[0] = v[0] - qe;
        max_H = H[0];
        max_t = 0;
      }
      // update ez
      if (en == tlen - 1 && H[en] > ez->mte)
        ez->mte = H[en], ez->mte_q = r - en;
      if (r - st == qlen - 1 && H[st] > ez->mqe)
        ez->mqe = H[st], ez->mqe_t = st;
      if (my_ksw_apply_zdrop(ez, 1, max_H, r, max_t, zdrop, e2, idx)) {
        //printf("zdrop in ksw_extd\n"); 
        break;
      }
      if (r == qlen + tlen - 2 && en == tlen - 1)
        ez->score = H[tlen - 1];
    } else { // find approximate max; Z-drop might be inaccurate, too.
      if (r > 0) {
        if (last_H0_t >= st && last_H0_t <= en && last_H0_t + 1 >= st && last_H0_t + 1 <= en) {
          int32_t d0 = v[last_H0_t];
          int32_t d1 = u[last_H0_t + 1];
          if (d0 > d1) H0 += d0;
          else H0 += d1, ++last_H0_t;
        } else if (last_H0_t >= st && last_H0_t <= en) {
          H0 += v[last_H0_t];
        } else {
          ++last_H0_t, H0 += u[last_H0_t];
        }
      } else H0 = v[0] - qe, last_H0_t = 0;
      if ((flag & KSW_EZ_APPROX_DROP) && my_ksw_apply_zdrop(ez, 1, H0, r, last_H0_t, zdrop, e2, idx)) break;
      if (r == qlen + tlen - 2 && en == tlen - 1)
        ez->score = H0;
    }
    last_st = st, last_en = en;
  }
  kfree(km, mem);
  if (!approx_max) kfree(km, H);
  if (with_cigar) { // backtrack
    int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
    if (!ez->zdropped && !(flag&KSW_EZ_EXTZ_ONLY)) {
      my_ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
    } else if (!ez->zdropped && (flag&KSW_EZ_EXTZ_ONLY) && ez->mqe + end_bonus > (int)ez->max) {
      ez->reach_end = 1;
      my_ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col, ez->mqe_t, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
    } else if (ez->max_t >= 0 && ez->max_q >= 0) {
      my_ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
    }
    kfree(km, mem2); kfree(km, off); kfree(km, p);
  }
}

void my_ksw_extd2_sse(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
           int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez, int idx)
{
#define __dp_code_block1 \
  z = _mm_load_si128(&s[t]); \
  xt1 = _mm_load_si128(&x[t]);                     /* xt1 <- x[r-1][t..t+15] */ \
  tmp = _mm_srli_si128(xt1, 15);                   /* tmp <- x[r-1][t+15] */ \
  xt1 = _mm_or_si128(_mm_slli_si128(xt1, 1), x1_); /* xt1 <- x[r-1][t-1..t+14] */ \
  x1_ = tmp; \
  vt1 = _mm_load_si128(&v[t]);                     /* vt1 <- v[r-1][t..t+15] */ \
  tmp = _mm_srli_si128(vt1, 15);                   /* tmp <- v[r-1][t+15] */ \
  vt1 = _mm_or_si128(_mm_slli_si128(vt1, 1), v1_); /* vt1 <- v[r-1][t-1..t+14] */ \
  v1_ = tmp; \
  a = _mm_add_epi8(xt1, vt1);                      /* a <- x[r-1][t-1..t+14] + v[r-1][t-1..t+14] */ \
  ut = _mm_load_si128(&u[t]);                      /* ut <- u[t..t+15] */ \
  b = _mm_add_epi8(_mm_load_si128(&y[t]), ut);     /* b <- y[r-1][t..t+15] + u[r-1][t..t+15] */ \
  x2t1= _mm_load_si128(&x2[t]); \
  tmp = _mm_srli_si128(x2t1, 15); \
  x2t1= _mm_or_si128(_mm_slli_si128(x2t1, 1), x21_); \
  x21_= tmp; \
  a2= _mm_add_epi8(x2t1, vt1); \
  b2= _mm_add_epi8(_mm_load_si128(&y2[t]), ut);

#define __dp_code_block2 \
  _mm_store_si128(&u[t], _mm_sub_epi8(z, vt1));    /* u[r][t..t+15] <- z - v[r-1][t-1..t+14] */ \
  _mm_store_si128(&v[t], _mm_sub_epi8(z, ut));     /* v[r][t..t+15] <- z - u[r-1][t..t+15] */ \
  tmp = _mm_sub_epi8(z, q_); \
  a = _mm_sub_epi8(a, tmp); \
  b = _mm_sub_epi8(b, tmp); \
  tmp = _mm_sub_epi8(z, q2_); \
  a2= _mm_sub_epi8(a2, tmp); \
  b2= _mm_sub_epi8(b2, tmp);

  int r, t, qe = q + e, n_col_, *off = 0, *off_end = 0, tlen_, qlen_, last_st, last_en, wl, wr, max_sc, min_sc, long_thres, long_diff;
  int with_cigar = !(flag&KSW_EZ_SCORE_ONLY), approx_max = !!(flag&KSW_EZ_APPROX_MAX);
  int32_t *H = 0, H0 = 0, last_H0_t = 0;
  uint8_t *qr, *sf, *mem, *mem2 = 0;
  __m128i q_, q2_, qe_, qe2_, zero_, sc_mch_, sc_mis_, m1_, sc_N_;
  __m128i *u, *v, *x, *y, *x2, *y2, *s, *p = 0;

  my_ksw_reset_extz(ez);
  if (m <= 1 || qlen <= 0 || tlen <= 0) return;

  if (q2 + e2 < q + e) t = q, q = q2, q2 = t, t = e, e = e2, e2 = t; // make sure q+e no larger than q2+e2

  zero_   = _mm_set1_epi8(0);
  q_      = _mm_set1_epi8(q);
  q2_     = _mm_set1_epi8(q2);
  qe_     = _mm_set1_epi8(q + e);
  qe2_    = _mm_set1_epi8(q2 + e2);
  sc_mch_ = _mm_set1_epi8(mat[0]);
  sc_mis_ = _mm_set1_epi8(mat[1]);
  sc_N_   = mat[m*m-1] == 0? _mm_set1_epi8(-e2) : _mm_set1_epi8(mat[m*m-1]);
  m1_     = _mm_set1_epi8(m - 1); // wildcard

  if (w < 0) w = tlen > qlen? tlen : qlen;
  wl = wr = w;
  tlen_ = (tlen + 15) / 16;
  n_col_ = qlen < tlen? qlen : tlen;
  n_col_ = ((n_col_ < w + 1? n_col_ : w + 1) + 15) / 16 + 1;
  qlen_ = (qlen + 15) / 16;
  for (t = 1, max_sc = mat[0], min_sc = mat[1]; t < m * m; ++t) {
    max_sc = max_sc > mat[t]? max_sc : mat[t];
    min_sc = min_sc < mat[t]? min_sc : mat[t];
  }
  if (-min_sc > 2 * (q + e)) return; // otherwise, we won't see any mismatches

  long_thres = e != e2? (q2 - q) / (e - e2) - 1 : 0;
  if (q2 + e2 + long_thres * e2 > q + e + long_thres * e)
    ++long_thres;
  long_diff = long_thres * (e - e2) - (q2 - q) - e2;

  mem = (uint8_t*)kcalloc(km, tlen_ * 8 + qlen_ + 1, 16);
  u = (__m128i*)(((size_t)mem + 15) >> 4 << 4); // 16-byte aligned
  v = u + tlen_, x = v + tlen_, y = x + tlen_, x2 = y + tlen_, y2 = x2 + tlen_;
  s = y2 + tlen_, sf = (uint8_t*)(s + tlen_), qr = sf + tlen_ * 16;
  memset(u,  -q  - e,  tlen_ * 16);
  memset(v,  -q  - e,  tlen_ * 16);
  memset(x,  -q  - e,  tlen_ * 16);
  memset(y,  -q  - e,  tlen_ * 16);
  memset(x2, -q2 - e2, tlen_ * 16);
  memset(y2, -q2 - e2, tlen_ * 16);
  if (!approx_max) {
    H = (int32_t*)kmalloc(km, tlen_ * 16 * 4);
    for (t = 0; t < tlen_ * 16; ++t) H[t] = KSW_NEG_INF;
  }
  if (with_cigar) {
    mem2 = (uint8_t*)kmalloc(km, ((size_t)(qlen + tlen - 1) * n_col_ + 1) * 16);
    p = (__m128i*)(((size_t)mem2 + 15) >> 4 << 4);
    off = (int*)kmalloc(km, (qlen + tlen - 1) * sizeof(int) * 2);
    off_end = off + qlen + tlen - 1;
  }

  for (t = 0; t < qlen; ++t) qr[t] = query[qlen - 1 - t];
  memcpy(sf, target, tlen);

  for (r = 0, last_st = last_en = -1; r < qlen + tlen - 1; ++r) {
    int st = 0, en = tlen - 1, st0, en0, st_, en_;
    int8_t x1, x21, v1;
    uint8_t *qrr = qr + (qlen - 1 - r);
    int8_t *u8 = (int8_t*)u, *v8 = (int8_t*)v, *x8 = (int8_t*)x, *x28 = (int8_t*)x2;
    __m128i x1_, x21_, v1_;
    // find the boundaries
    if (st < r - qlen + 1) st = r - qlen + 1;
    if (en > r) en = r;
    if (st < (r-wr+1)>>1) st = (r-wr+1)>>1; // take the ceil
    if (en > (r+wl)>>1) en = (r+wl)>>1; // take the floor
    if (st > en) {
      printf("zdropped in st > en\n");
      ez->zdropped = 1;
      break;
    }
    st0 = st, en0 = en;
    st = st / 16 * 16, en = (en + 16) / 16 * 16 - 1;
    // set boundary conditions
    if (st > 0) {
      if (st - 1 >= last_st && st - 1 <= last_en) {
        x1 = x8[st - 1], x21 = x28[st - 1], v1 = v8[st - 1]; // (r-1,s-1) calculated in the last round
      } else {
        x1 = -q - e, x21 = -q2 - e2;
        v1 = -q - e;
      }
    } else {
      x1 = -q - e, x21 = -q2 - e2;
      v1 = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
    }
    if (en >= r) {
      ((int8_t*)y)[r] = -q - e, ((int8_t*)y2)[r] = -q2 - e2;
      u8[r] = r == 0? -q - e : r < long_thres? -e : r == long_thres? long_diff : -e2;
    }
    // loop fission: set scores first
    if (!(flag & KSW_EZ_GENERIC_SC)) {
      for (t = st0; t <= en0; t += 16) {
        __m128i sq, st, tmp, mask;
        sq = _mm_loadu_si128((__m128i*)&sf[t]);
        st = _mm_loadu_si128((__m128i*)&qrr[t]);
        mask = _mm_or_si128(_mm_cmpeq_epi8(sq, m1_), _mm_cmpeq_epi8(st, m1_));
        tmp = _mm_cmpeq_epi8(sq, st);
#ifdef __SSE4_1__
        tmp = _mm_blendv_epi8(sc_mis_, sc_mch_, tmp);
        tmp = _mm_blendv_epi8(tmp,     sc_N_,   mask);
#else
        tmp = _mm_or_si128(_mm_andnot_si128(tmp,  sc_mis_), _mm_and_si128(tmp,  sc_mch_));
        tmp = _mm_or_si128(_mm_andnot_si128(mask, tmp),     _mm_and_si128(mask, sc_N_));
#endif
        _mm_storeu_si128((__m128i*)((int8_t*)s + t), tmp);
      }
    } else {
      for (t = st0; t <= en0; ++t)
        ((uint8_t*)s)[t] = mat[sf[t] * m + qrr[t]];
    }
    // core loop
    x1_  = _mm_cvtsi32_si128((uint8_t)x1);
    x21_ = _mm_cvtsi32_si128((uint8_t)x21);
    v1_  = _mm_cvtsi32_si128((uint8_t)v1);
    st_ = st / 16, en_ = en / 16;
    assert(en_ - st_ + 1 <= n_col_);
    if (!with_cigar) { // score only
      for (t = st_; t <= en_; ++t) {
        __m128i z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
        __dp_code_block1;
#ifdef __SSE4_1__
        z = _mm_max_epi8(z, a);
        z = _mm_max_epi8(z, b);
        z = _mm_max_epi8(z, a2);
        z = _mm_max_epi8(z, b2);
        z = _mm_min_epi8(z, sc_mch_);
        __dp_code_block2; // save u[] and v[]; update a, b, a2 and b2
        _mm_store_si128(&x[t],  _mm_sub_epi8(_mm_max_epi8(a,  zero_), qe_));
        _mm_store_si128(&y[t],  _mm_sub_epi8(_mm_max_epi8(b,  zero_), qe_));
        _mm_store_si128(&x2[t], _mm_sub_epi8(_mm_max_epi8(a2, zero_), qe2_));
        _mm_store_si128(&y2[t], _mm_sub_epi8(_mm_max_epi8(b2, zero_), qe2_));
#else
        tmp = _mm_cmpgt_epi8(a,  z);
        z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a));
        tmp = _mm_cmpgt_epi8(b,  z);
        z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b));
        tmp = _mm_cmpgt_epi8(a2, z);
        z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a2));
        tmp = _mm_cmpgt_epi8(b2, z);
        z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b2));
        tmp = _mm_cmplt_epi8(sc_mch_, z);
        z = _mm_or_si128(_mm_and_si128(tmp, sc_mch_), _mm_andnot_si128(tmp, z));
        __dp_code_block2;
        tmp = _mm_cmpgt_epi8(a, zero_);
        _mm_store_si128(&x[t],  _mm_sub_epi8(_mm_and_si128(tmp, a),  qe_));
        tmp = _mm_cmpgt_epi8(b, zero_);
        _mm_store_si128(&y[t],  _mm_sub_epi8(_mm_and_si128(tmp, b),  qe_));
        tmp = _mm_cmpgt_epi8(a2, zero_);
        _mm_store_si128(&x2[t], _mm_sub_epi8(_mm_and_si128(tmp, a2), qe2_));
        tmp = _mm_cmpgt_epi8(b2, zero_);
        _mm_store_si128(&y2[t], _mm_sub_epi8(_mm_and_si128(tmp, b2), qe2_));
#endif
      }
    } else if (!(flag&KSW_EZ_RIGHT)) { // gap left-alignment
      __m128i *pr = p + (size_t)r * n_col_ - st_;
      off[r] = st, off_end[r] = en;
      for (t = st_; t <= en_; ++t) {
        __m128i d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
        __dp_code_block1;
#ifdef __SSE4_1__
        d = _mm_and_si128(_mm_cmpgt_epi8(a, z), _mm_set1_epi8(1));       // d = a  > z? 1 : 0
        z = _mm_max_epi8(z, a);
        d = _mm_blendv_epi8(d, _mm_set1_epi8(2), _mm_cmpgt_epi8(b,  z)); // d = b  > z? 2 : d
        z = _mm_max_epi8(z, b);
        d = _mm_blendv_epi8(d, _mm_set1_epi8(3), _mm_cmpgt_epi8(a2, z)); // d = a2 > z? 3 : d
        z = _mm_max_epi8(z, a2);
        d = _mm_blendv_epi8(d, _mm_set1_epi8(4), _mm_cmpgt_epi8(b2, z)); // d = a2 > z? 3 : d
        z = _mm_max_epi8(z, b2);
        z = _mm_min_epi8(z, sc_mch_);
#else // we need to emulate SSE4.1 intrinsics _mm_max_epi8() and _mm_blendv_epi8()
        tmp = _mm_cmpgt_epi8(a,  z);
        d = _mm_and_si128(tmp, _mm_set1_epi8(1));
        z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a));
        tmp = _mm_cmpgt_epi8(b,  z);
        d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, _mm_set1_epi8(2)));
        z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b));
        tmp = _mm_cmpgt_epi8(a2, z);
        d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, _mm_set1_epi8(3)));
        z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, a2));
        tmp = _mm_cmpgt_epi8(b2, z);
        d = _mm_or_si128(_mm_andnot_si128(tmp, d), _mm_and_si128(tmp, _mm_set1_epi8(4)));
        z = _mm_or_si128(_mm_andnot_si128(tmp, z), _mm_and_si128(tmp, b2));
        tmp = _mm_cmplt_epi8(sc_mch_, z);
        z = _mm_or_si128(_mm_and_si128(tmp, sc_mch_), _mm_andnot_si128(tmp, z));
#endif
        __dp_code_block2;
        tmp = _mm_cmpgt_epi8(a, zero_);
        _mm_store_si128(&x[t],  _mm_sub_epi8(_mm_and_si128(tmp, a),  qe_));
        d = _mm_or_si128(d, _mm_and_si128(tmp, _mm_set1_epi8(0x08))); // d = a > 0? 1<<3 : 0
        tmp = _mm_cmpgt_epi8(b, zero_);
        _mm_store_si128(&y[t],  _mm_sub_epi8(_mm_and_si128(tmp, b),  qe_));
        d = _mm_or_si128(d, _mm_and_si128(tmp, _mm_set1_epi8(0x10))); // d = b > 0? 1<<4 : 0
        tmp = _mm_cmpgt_epi8(a2, zero_);
        _mm_store_si128(&x2[t], _mm_sub_epi8(_mm_and_si128(tmp, a2), qe2_));
        d = _mm_or_si128(d, _mm_and_si128(tmp, _mm_set1_epi8(0x20))); // d = a > 0? 1<<5 : 0
        tmp = _mm_cmpgt_epi8(b2, zero_);
        _mm_store_si128(&y2[t], _mm_sub_epi8(_mm_and_si128(tmp, b2), qe2_));
        d = _mm_or_si128(d, _mm_and_si128(tmp, _mm_set1_epi8(0x40))); // d = b > 0? 1<<6 : 0
        _mm_store_si128(&pr[t], d);
      }
    } else { // gap right-alignment
      __m128i *pr = p + (size_t)r * n_col_ - st_;
      off[r] = st, off_end[r] = en;
      for (t = st_; t <= en_; ++t) {
        __m128i d, z, a, b, a2, b2, xt1, x2t1, vt1, ut, tmp;
        __dp_code_block1;
#ifdef __SSE4_1__
        d = _mm_andnot_si128(_mm_cmpgt_epi8(z, a), _mm_set1_epi8(1));    // d = z > a?  0 : 1
        z = _mm_max_epi8(z, a);
        d = _mm_blendv_epi8(_mm_set1_epi8(2), d, _mm_cmpgt_epi8(z, b));  // d = z > b?  d : 2
        z = _mm_max_epi8(z, b);
        d = _mm_blendv_epi8(_mm_set1_epi8(3), d, _mm_cmpgt_epi8(z, a2)); // d = z > a2? d : 3
        z = _mm_max_epi8(z, a2);
        d = _mm_blendv_epi8(_mm_set1_epi8(4), d, _mm_cmpgt_epi8(z, b2)); // d = z > b2? d : 4
        z = _mm_max_epi8(z, b2);
        z = _mm_min_epi8(z, sc_mch_);
#else // we need to emulate SSE4.1 intrinsics _mm_max_epi8() and _mm_blendv_epi8()
        tmp = _mm_cmpgt_epi8(z, a);
        d = _mm_andnot_si128(tmp, _mm_set1_epi8(1));
        z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, a));
        tmp = _mm_cmpgt_epi8(z, b);
        d = _mm_or_si128(_mm_and_si128(tmp, d), _mm_andnot_si128(tmp, _mm_set1_epi8(2)));
        z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, b));
        tmp = _mm_cmpgt_epi8(z, a2);
        d = _mm_or_si128(_mm_and_si128(tmp, d), _mm_andnot_si128(tmp, _mm_set1_epi8(3)));
        z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, a2));
        tmp = _mm_cmpgt_epi8(z, b2);
        d = _mm_or_si128(_mm_and_si128(tmp, d), _mm_andnot_si128(tmp, _mm_set1_epi8(4)));
        z = _mm_or_si128(_mm_and_si128(tmp, z), _mm_andnot_si128(tmp, b2));
        tmp = _mm_cmplt_epi8(sc_mch_, z);
        z = _mm_or_si128(_mm_and_si128(tmp, sc_mch_), _mm_andnot_si128(tmp, z));
#endif
        __dp_code_block2;
        tmp = _mm_cmpgt_epi8(zero_, a);
        _mm_store_si128(&x[t],  _mm_sub_epi8(_mm_andnot_si128(tmp, a),  qe_));
        d = _mm_or_si128(d, _mm_andnot_si128(tmp, _mm_set1_epi8(0x08))); // d = a > 0? 1<<3 : 0
        tmp = _mm_cmpgt_epi8(zero_, b);
        _mm_store_si128(&y[t],  _mm_sub_epi8(_mm_andnot_si128(tmp, b),  qe_));
        d = _mm_or_si128(d, _mm_andnot_si128(tmp, _mm_set1_epi8(0x10))); // d = b > 0? 1<<4 : 0
        tmp = _mm_cmpgt_epi8(zero_, a2);
        _mm_store_si128(&x2[t], _mm_sub_epi8(_mm_andnot_si128(tmp, a2), qe2_));
        d = _mm_or_si128(d, _mm_andnot_si128(tmp, _mm_set1_epi8(0x20))); // d = a > 0? 1<<5 : 0
        tmp = _mm_cmpgt_epi8(zero_, b2);
        _mm_store_si128(&y2[t], _mm_sub_epi8(_mm_andnot_si128(tmp, b2), qe2_));
        d = _mm_or_si128(d, _mm_andnot_si128(tmp, _mm_set1_epi8(0x40))); // d = b > 0? 1<<6 : 0
        _mm_store_si128(&pr[t], d);
      }
    }
    if (!approx_max) { // find the exact max with a 32-bit score array
      int32_t max_H, max_t;
      // compute H[], max_H and max_t
      if (r > 0) {
        int32_t HH[4], tt[4], en1 = st0 + (en0 - st0) / 4 * 4, i;
        __m128i max_H_, max_t_;
        max_H = H[en0] = en0 > 0? H[en0-1] + u8[en0] : H[en0] + v8[en0]; // special casing the last element
        max_t = en0;
        max_H_ = _mm_set1_epi32(max_H);
        max_t_ = _mm_set1_epi32(max_t);
        for (t = st0; t < en1; t += 4) { // this implements: H[t]+=v8[t]-qe; if(H[t]>max_H) max_H=H[t],max_t=t;
          __m128i H1, tmp, t_;
          H1 = _mm_loadu_si128((__m128i*)&H[t]);
          t_ = _mm_setr_epi32(v8[t], v8[t+1], v8[t+2], v8[t+3]);
          H1 = _mm_add_epi32(H1, t_);
          _mm_storeu_si128((__m128i*)&H[t], H1);
          t_ = _mm_set1_epi32(t);
          tmp = _mm_cmpgt_epi32(H1, max_H_);
#ifdef __SSE4_1__
          max_H_ = _mm_blendv_epi8(max_H_, H1, tmp);
          max_t_ = _mm_blendv_epi8(max_t_, t_, tmp);
#else
          max_H_ = _mm_or_si128(_mm_and_si128(tmp, H1), _mm_andnot_si128(tmp, max_H_));
          max_t_ = _mm_or_si128(_mm_and_si128(tmp, t_), _mm_andnot_si128(tmp, max_t_));
#endif
        }
        _mm_storeu_si128((__m128i*)HH, max_H_);
        _mm_storeu_si128((__m128i*)tt, max_t_);
        for (i = 0; i < 4; ++i)
          if (max_H < HH[i]) max_H = HH[i], max_t = tt[i] + i;
        for (; t < en0; ++t) { // for the rest of values that haven't been computed with SSE
          H[t] += (int32_t)v8[t];
          if (H[t] > max_H)
            max_H = H[t], max_t = t;
        }
      } else H[0] = v8[0] - qe, max_H = H[0], max_t = 0; // special casing r==0
      // update ez
      if (en0 == tlen - 1 && H[en0] > ez->mte)
        ez->mte = H[en0], ez->mte_q = r - en;
      if (r - st0 == qlen - 1 && H[st0] > ez->mqe)
        ez->mqe = H[st0], ez->mqe_t = st0;
      if (my_ksw_apply_zdrop(ez, 1, max_H, r, max_t, zdrop, e2, idx)) {
        printf("zdrop in ksw_extd\n"); 
        break;
      }
      if (r == qlen + tlen - 2 && en0 == tlen - 1)
        ez->score = H[tlen - 1];
    } else { // find approximate max; Z-drop might be inaccurate, too.
      if (r > 0) {
        if (last_H0_t >= st0 && last_H0_t <= en0 && last_H0_t + 1 >= st0 && last_H0_t + 1 <= en0) {
          int32_t d0 = v8[last_H0_t];
          int32_t d1 = u8[last_H0_t + 1];
          if (d0 > d1) H0 += d0;
          else H0 += d1, ++last_H0_t;
        } else if (last_H0_t >= st0 && last_H0_t <= en0) {
          H0 += v8[last_H0_t];
        } else {
          ++last_H0_t, H0 += u8[last_H0_t];
        }
      } else H0 = v8[0] - qe, last_H0_t = 0;
      if ((flag & KSW_EZ_APPROX_DROP) && my_ksw_apply_zdrop(ez, 1, H0, r, last_H0_t, zdrop, e2, idx)) break;
      if (r == qlen + tlen - 2 && en0 == tlen - 1)
        ez->score = H0;
    }
    last_st = st, last_en = en;
    //for (t = st0; t <= en0; ++t) printf("(%d,%d)\t(%d,%d,%d,%d)\t%d\n", r, t, ((int8_t*)u)[t], ((int8_t*)v)[t], ((int8_t*)x)[t], ((int8_t*)y)[t], H[t]); // for debugging
  }
  kfree(km, mem);
  if (!approx_max) kfree(km, H);
  if (with_cigar) { // backtrack
    int rev_cigar = !!(flag & KSW_EZ_REV_CIGAR);
    if (!ez->zdropped && !(flag&KSW_EZ_EXTZ_ONLY)) {
      my_ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, tlen-1, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
    } else if (!ez->zdropped && (flag&KSW_EZ_EXTZ_ONLY) && ez->mqe + end_bonus > (int)ez->max) {
      ez->reach_end = 1;
      my_ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, ez->mqe_t, qlen-1, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
    } else if (ez->max_t >= 0 && ez->max_q >= 0) {
      my_ksw_backtrack(km, 1, rev_cigar, 0, (uint8_t*)p, off, off_end, n_col_*16, ez->max_t, ez->max_q, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
    }
    kfree(km, mem2); kfree(km, off);
  }
}


static void mm_align_pair(void *km, const mm_mapopt_t *opt, int qlen, const uint8_t *qseq, int tlen, const uint8_t *tseq, const int8_t *mat, int w, int end_bonus, int zdrop, int flag, ksw_extz_t *ez, double& cur_cells, double& sw_time, int idx)
{

  struct timespec time1, time2, time_diff;
  cur_cells += (double) (qlen * tlen);
  clock_gettime(CLOCK_REALTIME, &time1);
  ksw_gg2(km, qlen, qseq, tlen, tseq, 5, mat, opt->q, opt->e, opt->q2, opt->e2, w, zdrop, end_bonus, flag, ez, idx);
  clock_gettime(CLOCK_REALTIME, &time2);
  // time_diff = diff_time(time1, time2);
  // sw_time += (long)(time_diff.tv_sec * 1e9 + time_diff.tv_nsec);
}

static void mm_append_cigar(mm_reg1_t *r, uint32_t n_cigar, uint32_t *cigar) // TODO: this calls the libc realloc()
{
  mm_extra_t *p;
  if (n_cigar == 0) return;
  if (r->p == 0) {
    uint32_t capacity = n_cigar + sizeof(mm_extra_t)/4;
    kroundup32(capacity);
    r->p = (mm_extra_t*)calloc(capacity, 4);
    r->p->capacity = capacity;
  } else if (r->p->n_cigar + n_cigar + sizeof(mm_extra_t)/4 > r->p->capacity) {
    r->p->capacity = r->p->n_cigar + n_cigar + sizeof(mm_extra_t)/4;
    kroundup32(r->p->capacity);
    r->p = (mm_extra_t*)realloc(r->p, r->p->capacity * 4);
  }
  p = r->p;
  if (p->n_cigar > 0 && (p->cigar[p->n_cigar-1]&0xf) == (cigar[0]&0xf)) { // same CIGAR op at the boundary
    p->cigar[p->n_cigar-1] += cigar[0]>>4<<4;
    if (n_cigar > 1) memcpy(p->cigar + p->n_cigar, cigar + 1, (n_cigar - 1) * 4);
    p->n_cigar += n_cigar - 1;
  } else {
    memcpy(p->cigar + p->n_cigar, cigar, n_cigar * 4);
    p->n_cigar += n_cigar;
  }
}

static void mm_update_extra(mm_reg1_t *r, const uint8_t *qseq, const uint8_t *tseq, const int8_t *mat, int8_t q, int8_t e)
{
  uint32_t k, l;
  int32_t s = 0, max = 0, qshift, tshift, toff = 0, qoff = 0;
  mm_extra_t *p = r->p;
  if (p == 0) return;
  mm_fix_cigar(r, qseq, tseq, &qshift, &tshift);
  qseq += qshift, tseq += tshift; // qseq and tseq may be shifted due to the removal of leading I/D
  r->blen = r->mlen = 0;
  for (k = 0; k < p->n_cigar; ++k) {
    uint32_t op = p->cigar[k]&0xf, len = p->cigar[k]>>4;
    if (op == 0) { // match/mismatch
      int n_ambi = 0, n_diff = 0;
      for (l = 0; l < len; ++l) {
        int cq = qseq[qoff + l], ct = tseq[toff + l];
        if (ct > 3 || cq > 3) ++n_ambi;
        else if (ct != cq) ++n_diff;
        s += mat[ct * 5 + cq];
        if (s < 0) s = 0;
        else max = max > s? max : s;
      }
      r->blen += len - n_ambi, r->mlen += len - (n_ambi + n_diff), p->n_ambi += n_ambi;
      toff += len, qoff += len;
    } else if (op == 1) { // insertion
      int n_ambi = 0;
      for (l = 0; l < len; ++l)
        if (qseq[qoff + l] > 3) ++n_ambi;
      r->blen += len - n_ambi, p->n_ambi += n_ambi;
      s -= q + e * len;
      if (s < 0) s = 0;
      qoff += len;
    } else if (op == 2) { // deletion
      int n_ambi = 0;
      for (l = 0; l < len; ++l)
        if (tseq[toff + l] > 3) ++n_ambi;
      r->blen += len - n_ambi, p->n_ambi += n_ambi;
      s -= q + e * len;
      if (s < 0) s = 0;
      toff += len;
    } else if (op == 3) { // intron
      toff += len;
    }
  }
  p->dp_max = max;
}

static inline void update_max_zdrop(int32_t score, int i, int j, int32_t *max, int *max_i, int *max_j, int e, int *max_zdrop, int pos[2][2])
{
  if (score < *max) {
    int li = i - *max_i;
    int lj = j - *max_j;
    int diff = li > lj? li - lj : lj - li;
    int z = *max - score - diff * e;
    if (z > *max_zdrop) {
      *max_zdrop = z;
      pos[0][0] = *max_i, pos[0][1] = i + 1;
      pos[1][0] = *max_j, pos[1][1] = j + 1;
    }
  } else *max = score, *max_i = i, *max_j = j;
}

static int mm_test_zdrop(void *km, const mm_mapopt_t *opt, const uint8_t *qseq, const uint8_t *tseq, uint32_t n_cigar, uint32_t *cigar, const int8_t *mat)
{
  uint32_t k;
  int32_t score = 0, max = INT32_MIN, max_i = -1, max_j = -1, i = 0, j = 0, max_zdrop = 0;
  int pos[2][2] = {{-1, -1}, {-1, -1}}, q_len, t_len;

  // find the score and the region where score drops most along diagonal
  for (k = 0, score = 0; k < n_cigar; ++k) {
    uint32_t l, op = cigar[k]&0xf, len = cigar[k]>>4;
    if (op == 0) {
      for (l = 0; l < len; ++l) {
        score += mat[tseq[i + l] * 5 + qseq[j + l]];
        update_max_zdrop(score, i+l, j+l, &max, &max_i, &max_j, opt->e, &max_zdrop, pos);
      }
      i += len, j += len;
    } else if (op == 1 || op == 2 || op == 3) {
      score -= opt->q + opt->e * len;
      if (op == 1) j += len; // insertion
      else i += len;         // deletion
      update_max_zdrop(score, i, j, &max, &max_i, &max_j, opt->e, &max_zdrop, pos);
    }
  }

  // test if there is an inversion in the most dropped region
  q_len = pos[1][1] - pos[1][0], t_len = pos[0][1] - pos[0][0];

  if (!(opt->flag&(MM_F_SPLICE|MM_F_SR|MM_F_FOR_ONLY|MM_F_REV_ONLY)) && max_zdrop > opt->zdrop_inv && q_len < opt->max_gap && t_len < opt->max_gap) {
    uint8_t *qseq2;
    void *qp;
    int q_off, t_off;
    qseq2 = (uint8_t*)kmalloc(km, q_len);
    for (i = 0; i < q_len; ++i) {
      int c = qseq[pos[1][1] - i - 1];
      qseq2[i] = c >= 4? 4 : 3 - c;
    }
    qp = ksw_ll_qinit(km, 2, q_len, qseq2, 5, mat);
    score = ksw_ll_i16(qp, t_len, tseq + pos[0][0], opt->q, opt->e, &q_off, &t_off);
    kfree(km, qseq2);
    kfree(km, qp);
    if (score >= opt->min_chain_score * opt->a && score >= opt->min_dp_max)
      return 2; // there is a potential inversion
  }
  return max_zdrop > opt->zdrop? 1 : 0;
}



void mm_align_baseline(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, uint8_t *qseq0[2], mm_reg1_t *r, mm_reg1_t *r2, \
    int n_a, mm128_t *a, ksw_extz_t *ez, int splice_flag, double& cur_cells, double& sw_time, int idx)
{
  int32_t rid = a[r->as].x<<1>>33, rev = a[r->as].x>>63, as1, cnt1;
  uint8_t *tseq, *qseq;
  int32_t i, l, bw, dropped = 0, extra_flag = 0, rs0, re0, qs0, qe0;
  int32_t rs, re, qs, qe;
  int32_t rs1, qs1, re1, qe1;
  int8_t mat[25];

  r2->cnt = 0;
  if (r->cnt == 0) return;
  ksw_gen_simple_mat(5, mat, opt->a, opt->b, opt->sc_ambi);
  bw = (int)(opt->bw * 1.5 + 1.);
  mm_max_stretch(r, a, &as1, &cnt1);
  rs = (int32_t)a[as1].x + 1 - (int32_t)(a[as1].y>>32&0xff);
  qs = (int32_t)a[as1].y + 1 - (int32_t)(a[as1].y>>32&0xff);
  re = (int32_t)a[as1+cnt1-1].x + 1;
  qe = (int32_t)a[as1+cnt1-1].y + 1;

  /* Look for the start and end of regions to perform DP. This sounds easy
   * but is in fact tricky. Excessively small regions lead to unnecessary
   * clippings and lose alignable sequences. Excessively large regions
   * occasionally lead to large overlaps between two chains and may cause
   * loss of alignments in corner cases. */
  qs0 = 0, qe0 = qlen;
  l = qs;
  l += l * opt->a + opt->end_bonus > opt->q? (l * opt->a + opt->end_bonus - opt->q) / opt->e : 0;
  l = l < opt->bw? l : opt->bw;
  rs0 = rs - l > 0? rs - l : 0;
  l = qlen - qe;
  l += l * opt->a + opt->end_bonus > opt->q? (l * opt->a + opt->end_bonus - opt->q) / opt->e : 0;
  l = l < opt->bw? l : opt->bw;
  re0 = re + l < (int32_t)mi->seq[rid].len? re + l : mi->seq[rid].len;
  if (a[r->as].y & MM_SEED_SELF) {
    int max_ext = r->qs > r->rs? r->qs - r->rs : r->rs - r->qs;
    if (r->rs - rs0 > max_ext) rs0 = r->rs - max_ext;
    if (r->qs - qs0 > max_ext) qs0 = r->qs - max_ext;
    max_ext = r->qe > r->re? r->qe - r->re : r->re - r->qe;
    if (re0 - r->re > max_ext) re0 = r->re + max_ext;
    if (qe0 - r->qe > max_ext) qe0 = r->qe + max_ext;
  }

  tseq = (uint8_t*)kmalloc(km, re0 - rs0);

  if (qs > 0 && rs > 0) { // left extension
    qseq = &qseq0[rev][qs0];
    my_mm_idx_getseq(mi, rid, rs0, rs, tseq);
    mm_seq_rev(qs - qs0, qseq);
    mm_seq_rev(rs - rs0, tseq);
    //printf(" %d left extension %d x %d\n", idx, qs - qs0, rs - rs0);
    mm_align_pair(km, opt, qs - qs0, qseq, rs - rs0, tseq, mat, bw, opt->end_bonus, r->split_inv? opt->zdrop_inv : opt->zdrop, extra_flag|KSW_EZ_EXTZ_ONLY|KSW_EZ_RIGHT|KSW_EZ_REV_CIGAR, ez, cur_cells, sw_time, idx);
    if (ez->n_cigar > 0) {
      mm_append_cigar(r, ez->n_cigar, ez->cigar);
      r->p->dp_score += ez->max;
    }
    rs1 = rs - (ez->reach_end? ez->mqe_t + 1 : ez->max_t + 1);
    qs1 = qs - (ez->reach_end? qs - qs0 : ez->max_q + 1);
    mm_seq_rev(qs - qs0, qseq);
  } else rs1 = rs, qs1 = qs;
  re1 = rs, qe1 = qs;

  i = cnt1 - 1;
  bool is_break = false;
  //for (i = cnt1 - 1; i < cnt1; ++i) { // gap filling
    //if ((a[as1+i].y & (MM_SEED_IGNORE|MM_SEED_TANDEM)) && i != cnt1 - 1) continue;
    re = (int32_t)a[as1 + i].x + 1;
    qe = (int32_t)a[as1 + i].y + 1;
    re1 = re, qe1 = qe;
    if (i == cnt1 - 1 || (a[as1+i].y&MM_SEED_LONG_JOIN) || (qe - qs >= opt->min_ksw_len && re - rs >= opt->min_ksw_len)) {
      int j, bw1 = bw, zdrop_code;
      if (a[as1+i].y & MM_SEED_LONG_JOIN)
        bw1 = qe - qs > re - rs? qe - qs : re - rs;
      // perform alignment
      qseq = &qseq0[rev][qs];
      my_mm_idx_getseq(mi, rid, rs, re, tseq);
      my_ksw_reset_extz(ez);
      for (j = 0, ez->score = 0; j < qe - qs; ++j) {
        if (qseq[j] >= 4 || tseq[j] >= 4) ez->score += opt->e2;
        else ez->score += qseq[j] == tseq[j]? opt->a : -opt->b;
      }
      ez->cigar = my_push_cigar(km, &ez->n_cigar, &ez->m_cigar, ez->cigar, 0, qe - qs);
      // test Z-drop and inversion Z-drop
      if ((zdrop_code = mm_test_zdrop(km, opt, qseq, tseq, ez->n_cigar, ez->cigar, mat)) != 0) {
        //printf("%d, gap filling %d x %d\n",idx, qe - qs, re - rs);
        mm_align_pair(km, opt, qe - qs, qseq, re - rs, tseq, mat, bw1, -1, zdrop_code == 2? opt->zdrop_inv : opt->zdrop, extra_flag, ez, cur_cells, sw_time, idx); // second pass: lift approximate
      }
      // update CIGAR
      if (ez->n_cigar > 0)
        mm_append_cigar(r, ez->n_cigar, ez->cigar);
      if (ez->zdropped) { // truncated by Z-drop; TODO: sometimes Z-drop kicks in because the next seed placement is wrong. This can be fixed in principle.
        for (j = i - 1; j >= 0; --j)
          if ((int32_t)a[as1 + j].x <= rs + ez->max_t) {
            is_break = true;
            //break;
          }
        dropped = 1;
        if (j < 0) j = 0;
        r->p->dp_score += ez->max;
        re1 = rs + (ez->max_t + 1);
        qe1 = qs + (ez->max_q + 1);
        if (cnt1 - (j + 1) >= opt->min_cnt) {
          my_mm_split_reg(r, r2, as1 + j + 1 - r->as, qlen, a);
          if (zdrop_code == 2) r2->split_inv = 1;
        }
        is_break = true;
        //break;
      } else r->p->dp_score += ez->score;
      if (!is_break) {
        rs = re, qs = qe;
      }
    }
  //}

  if (!dropped && qe < qe0 && re < re0) { // right extension
    qseq = &qseq0[rev][qe];
    my_mm_idx_getseq(mi, rid, re, re0, tseq);
    //printf("%d, right extension %d x %d\n", idx, qe0 - qe, re0 - re);
    mm_align_pair(km, opt, qe0 - qe, qseq, re0 - re, tseq, mat, bw, opt->end_bonus, opt->zdrop, extra_flag|KSW_EZ_EXTZ_ONLY, ez, cur_cells, sw_time, idx);
    if (ez->n_cigar > 0) {
      mm_append_cigar(r, ez->n_cigar, ez->cigar);
      r->p->dp_score += ez->max;
    }
    re1 = re + (ez->reach_end? ez->mqe_t + 1 : ez->max_t + 1);
    qe1 = qe + (ez->reach_end? qe0 - qe : ez->max_q + 1);
  }
  //assert(qe1 <= qlen);

  r->rs = rs1, r->re = re1;
  if (rev) r->qs = qlen - qe1, r->qe = qlen - qs1;
  else r->qs = qs1, r->qe = qe1;

  //assert(re1 - rs1 <= re0 - rs0);
  if (r->p) {
    my_mm_idx_getseq(mi, rid, rs1, re1, tseq);
    mm_update_extra(r, &qseq0[r->rev][qs1], tseq, mat, opt->q, opt->e);
    if (rev && r->p->trans_strand)
      r->p->trans_strand ^= 3; // flip to the read strand
  }

  kfree(km, tseq);
}


void mmsw_baseline_compute(const mm_mapopt_t* opt, const mm_idx_t* mi, vector<align_input>& inputs, vector<align_output>& outputs, double& cur_cells, double& sw_time) {
  // void *km = alloc_km();
  uint32_t cigar[MAX_SEQ_LENGTH];
  for (int i = 0; i < inputs.size(); ++i) {
    align_output cur_output;
    ksw_extz_t ez;
     memset(&ez, 0, sizeof(ksw_extz_t));
    uint8_t *qseq_ptr[2];
    qseq_ptr[0] = inputs[i].qseq[0];
    qseq_ptr[1] = inputs[i].qseq[1];
    ez.cigar = cigar;
    ez.m_cigar = MAX_SEQ_LENGTH; 
    mm_align_baseline(NULL, opt, mi, inputs[i].qlen, qseq_ptr, \
        &(inputs[i].region[0].orig), &(inputs[i].region[1].orig), \
        inputs[i].n_a, inputs[i].a, &ez, 0, cur_cells, sw_time, i);
    cur_output.region[0] = inputs[i].region[0];   
    cur_output.region[1] = inputs[i].region[1];   
    cur_output.p = *(inputs[i].region[0].orig.p);
    assert(inputs[i].region[0].orig.p->n_cigar <= MAX_SEQ_LENGTH);
    memcpy(cur_output.cigar, inputs[i].region[0].orig.p->cigar, sizeof(uint32_t) * inputs[i].region[0].orig.p->n_cigar); 
    outputs.push_back(cur_output);   
  }
  // destroy_km(km);
}

