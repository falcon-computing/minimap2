#include <iostream>
#include <string>

#ifdef NDEBUG
#define LOG_HEADER "falcon-minimap2"
#endif
#include <glog/logging.h>
#include <gflags/gflags.h>
#include <gtest/gtest.h>

#include "minimap.h"
#include "mmpriv.h"
#include "MnmpGlobal.h"

#include "MnmpUtils.h"
#include "MnmpOptions.h"
#include "MnmpCpuStages.h"

#include "TestsCommon.h"

mm_mapopt_t *g_mnmpOpt;
mm_idxopt_t *g_mnmpIpt;
mm_idx_t *g_minimizer;
mm_idx_reader_t *g_idxReader;

extern int    g_numFp;
extern char **g_fn;

int    g_numFp;
char **g_fn;


int main(int argc, char *argv[]) {
  // Store original arguments
  std::stringstream l_cmdStr;
  for (int i = 0; i < argc; i++) {
    l_cmdStr << argv[i] << " ";
  }

  // Initialize Google Flags
  std::string l_version = "Falcon MiniMap2 Version: " + std::string(VERSION);
  google::SetVersionString(l_version.c_str());
  google::SetUsageMessage(argv[0]);
  google::ParseCommandLineFlags(&argc, &argv, true);

  // Initialize Google Log
  google::InitGoogleLogging(argv[0]);

  // Initialize Google Test
  ::testing::InitGoogleTest(&argc, argv);

  // Check arguments for input
  if (argc < 2) {
    LOG(ERROR) << "Unable to find input index file";
    return 1;
  }

  // Parse arguments for config
  g_mnmpOpt = (mm_mapopt_t*)malloc(sizeof(mm_mapopt_t));
  g_mnmpIpt = (mm_idxopt_t*)malloc(sizeof(mm_idxopt_t));

  int l_err = fc_set_opt();
  if (l_err != 0) {
    free(g_mnmpOpt);
    free(g_mnmpIpt);
    DLOG(ERROR) << "Unexpected arguments";
    return -1;
  }


  // Open index file
  std::string l_indexFile(argv[1]);
  std::string l_indexDumpFile(FLAGS_d);
  g_idxReader = mm_idx_reader_open(l_indexFile.c_str(), g_mnmpIpt, FLAGS_d.c_str());
  if (g_idxReader == 0) {
    mm_idx_reader_close(g_idxReader);
    g_idxReader = NULL;
    LOG(ERROR) << "Failed to open file " << l_indexFile;
    return 1;
  }
  if (!g_idxReader->is_idx  &&
      l_indexDumpFile == "" &&
      argc - 1 < 2            ) {
    mm_idx_reader_close(g_idxReader);
    g_idxReader = NULL;
    LOG(ERROR) << "Missing input: please specify a query file to map or option -d to keep the index";
    return 1;
  }
  if (g_mnmpOpt->best_n == 0         &&
      (g_mnmpOpt->flag & MM_F_CIGAR)   ) {
    LOG_IF(WARNING, VLOG_IS_ON(2)) << "`-N 0' reduces alignment accuracy. Please use --secondary=no to suppress secondary alignments."; 
  }
  g_minimizer = mm_idx_reader_read(g_idxReader, std::min(FLAGS_t, 3));
  if ((g_mnmpOpt->flag & MM_F_CIGAR)    &&
      (g_minimizer->flag & MM_I_NO_SEQ)   ) {
    mm_idx_destroy(g_minimizer);
    g_minimizer = NULL;
    mm_idx_reader_close(g_idxReader);
    g_idxReader = NULL;
    LOG(ERROR) << "The prebuilt index doesn't contain sequences.";
    return 1;
  }
  if ((g_mnmpOpt->flag & MM_F_OUT_SAM) &&
      g_idxReader->n_parts == 1          ) {
    // skip header
  }
  LOG_IF(INFO, VLOG_IS_ON(3)) << "Loaded/built the index for " << g_minimizer->n_seq << " target sequence(s)";
  if (argc != 2) {
    mm_mapopt_update(g_mnmpOpt, g_minimizer);
  }
  if (FLAGS_v >= 3) {
    mm_idx_stat(g_minimizer);
  }
  g_numFp = argc - 2;
  g_fn = &argv[2];

  int ret = RUN_ALL_TESTS();


  // Close index file 
  mm_idx_destroy(g_minimizer);
  mm_idx_reader_close(g_idxReader);

  // Clean up
  free(g_mnmpOpt);
  free(g_mnmpIpt);

  return ret;
}
