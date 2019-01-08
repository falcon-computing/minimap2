#include "mmsw_orig.h"
#include <string.h>

#include "MnmpWrapper.h"

// extern "C" {
//     void mm_align1(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int qlen, uint8_t *qseq0[2], mm_reg1_t *r, mm_reg1_t *r2, int n_a, mm128_t *a, ksw_extz_t *ez, int splice_flag);
// }

void mmsw_orig_compute(const mm_mapopt_t* opt, const mm_idx_t* mi, std::vector<align_input>& inputs, std::vector<align_output>& outputs) {
    // void *km = alloc_km();
    outputs.clear();
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
        mm_align1(NULL, opt, mi, inputs[i].qlen, qseq_ptr, \
                &(inputs[i].region[0].orig), &(inputs[i].region[1].orig), \
                inputs[i].n_a, inputs[i].a, &ez, opt->flag);
        cur_output.region[0] = inputs[i].region[0];   
        cur_output.region[1] = inputs[i].region[1];   
        cur_output.p = *(inputs[i].region[0].orig.p);
        memcpy(cur_output.cigar, inputs[i].region[0].orig.p->cigar, sizeof(uint32_t) * inputs[i].region[0].orig.p->n_cigar);
        // fix memory leakage
        free(inputs[i].region[0].orig.p);
        outputs.push_back(cur_output);
    }
    // destroy_km(km);
}






