#include <glog/logging.h>
#include <unordered_set>
#include <vector>

#include "blaze/Client.h"
#include "ksight/tools.h"
#include "MnmpSWWorker.h"
#include "../MnmpAlignFPGAUtils/host_interface/mmsw_orig.h"
#include "../MnmpAlignFPGAUtils/host_interface/mmsw_fpga.h"

#include "MnmpGlobal.h"

using namespace blaze;

#define mm_seq4_get(s, i)    ((s)[(i)>>3] >> (((i)&7)<<2) & 0xf)

static int my_mm_idx_getseq(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq)
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

static void mm_fix_cigar(mm_reg1_t *r, mm_extra_t *p, const uint8_t *qseq, const uint8_t *tseq, int *qshift, int *tshift)
{
	//mm_extra_t *p = r->p;
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

static void mm_update_extra(mm_reg1_t *r, mm_extra_t *p, const uint8_t *qseq, const uint8_t *tseq, const int8_t *mat, int8_t q, int8_t e)
{
	uint32_t k, l;
	int32_t s = 0, max = 0, qshift, tshift, toff = 0, qoff = 0;
	//mm_extra_t *p = r->p;
	if (p == 0) return;
	mm_fix_cigar(r, p, qseq, tseq, &qshift, &tshift);
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

MnmpSWWorker::MnmpSWWorker(
    MnmpSWClient* client,
    std::vector<align_input>& alignInputs,
    std::vector<align_output>& alignOutputs,
    int num_region): 
  client_(client),
  use_cpu_(false),
  input_(&alignInputs),
#ifndef DUPLICATE_OUTPUT
  output_(&alignOutputs),
#endif
  num_region_(num_region)
{
  PLACE_TIMER;

#ifdef DUPLICATE_OUTPUT
  output_.resize(num_region);
#else
  output_->resize(num_region);
#endif

  DLOG(INFO) << "fpga workload: num_region = " << num_region;
}

MnmpSWWorker:: ~MnmpSWWorker() {
  DLOG(INFO) << "Destroying worker";
}

// cpu compute routine
void MnmpSWWorker::compute() {
  ksight::AutoTimer __timer("compute on cpu");

  if (use_cpu_) {
#ifdef DUPLICATE_OUTPUT
    mmsw_orig_compute(g_mnmpOpt, g_minimizer, *input_, output_);
#else
    mmsw_orig_compute(g_mnmpOpt, g_minimizer, *input_, *output_);
#endif
  }
}

void MnmpSWWorker::getOutput(std::vector<align_output>& alignOutputs) {
  PLACE_TIMER;

  mm_mapopt_t *opt = g_mnmpOpt;
  mm_idx_t *mi = g_minimizer;;
  int8_t mat[25];
  ksw_gen_simple_mat(5, mat, opt->a, opt->b, opt->sc_ambi);
  for (int i = 0; i < num_region_; i++) {

    mm_reg1_t* r = &((*output_)[i].region[0].orig); 
    mm_reg1_t* r2 = &((*output_)[i].region[1].orig); 
    if ((*output_)[i].region[0].orig.split & 0x1) {
      r2->score = (int32_t)(r->score * ((float)r2->cnt / r->cnt) + .499);
      r->score -= r2->score;
    }
    mm_extra_t* p = &((*output_)[i].p);
    mm128_t* a = (*input_)[i].a;
    if (p->n_cigar > 0) {
      int32_t rs1 = r->rs;
      int32_t re1 = r->re;
      int32_t rev = a[r->as].x>>63;
      int32_t qs1 = rev ? (*input_)[i].qlen - r->qe : r->qs; 
      int32_t rid = a[r->as].x<<1>>33;
      uint8_t tseq[MAX_SEQ_LENGTH];
      my_mm_idx_getseq(mi, rid, rs1, re1, tseq);
      mm_update_extra(r, p, &(*input_)[i].qseq[r->rev][qs1], tseq, mat, opt->q, opt->e);
    }
  }
#ifdef DUPLICATE_OUTPUT
  // swap the buffers between worker and compute_stage
  std::swap(alignOutputs, output_);
#endif
} 

// currently all workload are calculated on fpga
void MnmpSWWorker::run() {
  DLOG(INFO) << "Computing mode: " << (use_cpu_ ? "CPU" : "CFX");
  if (use_cpu_) {
    compute();
  }
  else {
    //PLACE_TIMER1("compute on fpga");
    ksight::Timer timer;
    timer.start();

    // start a thread to run cpu
    // boost::thread t(boost::bind(&PairHMMWorker::compute, this));

    // set input
    {
      PLACE_TIMER1("serialize input");
#ifdef DUPLICATE_OUTPUT
      client_->setup(*input_, num_region_);
#else
      client_->setup(*input_, *output_, num_region_);
#endif
    }

    // start fpga run
    {
      PLACE_TIMER1("start()");
      client_->start();
    }

    // if cpu fallback is called, skip output copy
    {
      PLACE_TIMER1("copy output");

#ifdef DUPLICATE_OUTPUT
      // copy the data from MMSWTask output buffer
      fpga_outputs *results = (fpga_outputs *)client_->getOutputPtr(0); 
      memcpy(output_.data(), results, sizeof(align_output)*num_region_);
#else
      // get a pseudo-output token
      int ret = *((int *)client_->getOutputPtr(0));
      // assert(ret == 1);
#endif
    }
    // don't count cpu time
    ksight::ksight.add("compute on client", timer.stop());
    
    // t.join();
  }
}
