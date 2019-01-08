#ifndef MNMP_FLOW_GLOBAL_VARIABLES_H
#define MNMP_FLOW_GLOBAL_VARIABLES_H

#include <vector>

#include "minimap.h"
#include "htslib/sam.h"

extern mm_mapopt_t *g_mnmpOpt;
extern mm_idxopt_t *g_mnmpIpt;
extern mm_idx_t *g_minimizer;
extern mm_idx_reader_t *g_idxReader;
extern bam_hdr_t *g_bamHeader;

extern std::vector<mm_idx_t*> g_minimizerNumaList;

#endif
