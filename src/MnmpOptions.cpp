#include <gflags/gflags.h>
#include <boost/thread.hpp>
#include "MnmpOptions.h"
#include "MnmpGlobal.h"

#include "minimap.h"
#include "mmpriv.h"


DEFINE_string(d, "",
		"-d arg in original minimap2, dump index to file");
DEFINE_bool(a, true,
                "-a arg in original minimap2, output in SAM format");
DEFINE_bool(sam, true,
                "--sam arg in original minimap2, output in SAM format");
DEFINE_int32(t, 1,
                "-t arg in original minimap2, number of threads");
DEFINE_string(x, "sr",
                "-x arg in original minimap2, preset");
DEFINE_bool(frag, true,
                "--frag arg in original minimap2, enable frag g_mnmpOptde");
DEFINE_string(R, "",
                "-R arg in original minimap2, SAM read group line in a format like '@RG\\tID:foo\\tSM:bar'");

DEFINE_string(output_dir, "./",
                "output directory for SAM/BAM files");
DEFINE_int32(output_flag, -1,
                "Flag to specify output format: "
                "0: BAM (compressed); 1: BAM (uncompressed); 2: SAM");
DEFINE_int32(output_size, 128,
                "number of batches in a SAM/BAM file");
DEFINE_bool(inorder_output, false,
                "write all batches in order");
DEFINE_bool(sort, true,
                "apply coordinate sort");
DEFINE_bool(bam, false,
                "output in BAM format");
DEFINE_int32(extra_threads, 0,
                "number of extra threads to be reduces");

#ifdef BUILD_FPGA
DEFINE_bool(use_fpga, true,
                "enable fpga computation");
DEFINE_bool(fpga_only, false,
                "disable cpu for fpga-available stages");
DEFINE_int32(fpga_threads, 1,
                "number of FPGA threads");
DEFINE_string(fpga_path, "",
                "path to FPGA bitstream");
#ifdef LOCAL_BLAZE
DEFINE_string(blaze_conf, "",
                "config file for Blaze");
#endif
#endif

DEFINE_bool(use_numa, false,
                "use NUMA lib for hardware locality");

int fc_set_opt() {
  mm_idxopt_init(g_mnmpIpt);
  mm_mapopt_init(g_mnmpOpt);

  if (FLAGS_x == "sr" || FLAGS_x == "short") {
    g_mnmpIpt->flag = 0, g_mnmpIpt->k = 21, g_mnmpIpt->w = 11;
    g_mnmpOpt->flag |= MM_F_SR | MM_F_FRAG_MODE | MM_F_NO_PRINT_2ND | MM_F_2_IO_THREADS | MM_F_HEAP_SORT;
    g_mnmpOpt->pe_ori = 0<<1|1; // FR
    g_mnmpOpt->a = 2, g_mnmpOpt->b = 8, g_mnmpOpt->q = 12, g_mnmpOpt->e = 2, g_mnmpOpt->q2 = 24, g_mnmpOpt->e2 = 1;
    g_mnmpOpt->zdrop = g_mnmpOpt->zdrop_inv = 100;
    g_mnmpOpt->end_bonus = 10;
    g_mnmpOpt->max_frag_len = 800;
    g_mnmpOpt->max_gap = 100;
    g_mnmpOpt->bw = 100;
    g_mnmpOpt->pri_ratio = 0.5f;
    g_mnmpOpt->min_cnt = 2;
    g_mnmpOpt->min_chain_score = 25;
    g_mnmpOpt->min_dp_max = 40;
    g_mnmpOpt->best_n = 20;
    g_mnmpOpt->mid_occ = 1000;
    g_mnmpOpt->max_occ = 5000;
    g_mnmpOpt->mini_batch_size = 2560000;
  }
  else {
    return -1;
  }

  if (FLAGS_frag) {
    g_mnmpOpt->flag |= MM_F_FRAG_MODE;
  }
  else {
    g_mnmpOpt->flag &= ~MM_F_FRAG_MODE;
  }

  // Config output format
  if (FLAGS_a || FLAGS_sam) {
    g_mnmpOpt->flag |= MM_F_OUT_SAM | MM_F_CIGAR;
    FLAGS_output_flag = 2;
  }
  if (FLAGS_bam || FLAGS_output_flag == -1) { //Bam config will overwrite sam config
    FLAGS_output_flag = 1;
  } //TODO: add support for PAF format

  return 0;
}
