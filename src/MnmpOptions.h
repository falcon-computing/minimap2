#ifndef MNMP_FLOW_CONFIG_H
#define MNMP_FLOW_CONFIG_H

#include <string>
#include <gflags/gflags.h>

//extern mm_mapopt_t *g_mnmpOpt;
//extern mm_idxopt_t *g_mnmpIpt;


DECLARE_string(d);
DECLARE_bool(a);
DECLARE_bool(sam);
DECLARE_string(x);
DECLARE_bool(frag);
DECLARE_int32(t);
DECLARE_string(R);

DECLARE_string(output_dir);
DECLARE_int32(output_flag);
DECLARE_int32(output_size);
DECLARE_bool(inorder_output);
DECLARE_bool(sort);
DECLARE_bool(bam);

int fc_set_opt(void);

#endif
