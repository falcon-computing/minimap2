#include <iostream>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include <numa.h>

#ifdef NDEBUG
#define LOG_HEADER "falcon-minimap2"
#endif
#include <glog/logging.h>
#include <gflags/gflags.h>

#include "minimap.h"
#include "mmpriv.h"
#include "htslib/sam.h"
#include "MnmpGlobal.h"

#include "MnmpUtils.h"
#include "MnmpOptions.h"
#include "MnmpWrapper.h"
#include "MnmpCpuStages.h"

mm_mapopt_t *g_mnmpOpt;
mm_idxopt_t *g_mnmpIpt;
mm_idx_t *g_minimizer;
mm_idx_reader_t *g_idxReader;
bam_hdr_t *g_bamHeader;

std::vector<mm_idx_t*> g_minimizerNumaList;


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

  // Print arguments for records
  DLOG(INFO) << l_cmdStr.str();

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

  int l_numThreads = FLAGS_t;

  if (!boost::filesystem::create_directories(FLAGS_temp_dir)) {
    DLOG(INFO) << "Can not create temp folder: " << FLAGS_temp_dir;
  }

  // Start timer
  double l_startTime = realtime();

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

  int l_numNumaNodes = 1;
  if (numa_available() != -1) {
    l_numNumaNodes = numa_num_configured_nodes();
  }
  FLAGS_use_numa = (l_numNumaNodes > 1);
  DLOG_IF(INFO, FLAGS_use_numa) << "Found " << l_numNumaNodes << " NUMA nodes";

  struct bitmask *l_nodeMask = numa_parse_nodestring("0");
  if (FLAGS_use_numa)
    numa_set_membind(l_nodeMask);
  numa_free_nodemask(l_nodeMask);
  DLOG_IF(INFO, FLAGS_use_numa) << "Loading index on NUMA node 0";
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
  g_minimizerNumaList.push_back(g_minimizer);

  for (int l_nn = 1; l_nn < l_numNumaNodes; l_nn++) {
    l_nodeMask = numa_parse_nodestring(std::to_string(l_nn).c_str());
    numa_set_membind(l_nodeMask);
    numa_free_nodemask(l_nodeMask);
    DLOG_IF(INFO, FLAGS_use_numa) << "Loading index on NUMA node " << l_nn;
    mm_idx_reader_t *l_idxReader = mm_idx_reader_open(l_indexFile.c_str(), g_mnmpIpt, NULL);
    mm_idx_t *l_minimizer = mm_idx_reader_read(l_idxReader, std::min(FLAGS_t, 3));
    mm_idx_reader_close(l_idxReader);
    g_minimizerNumaList.push_back(l_minimizer);
  }
  if (FLAGS_use_numa)
    numa_set_membind(numa_all_nodes_ptr);

  LOG_IF(INFO, VLOG_IS_ON(3)) << "Loaded/built the index for " << g_minimizer->n_seq << " target sequence(s)";
  if (argc != 2) {
    mm_mapopt_update(g_mnmpOpt, g_minimizer);
  }
  else {
    return 0;
  }
  if (FLAGS_v >= 3) {
    mm_idx_stat(g_minimizer);
  }

  // Prepare header
  std::string l_headerStr = fc_write_sam_hdr(g_minimizer, FLAGS_R, VERSION, l_cmdStr.str());
  g_bamHeader = sam_hdr_parse(l_headerStr.length(), l_headerStr.c_str());
  g_bamHeader->l_text = l_headerStr.length();
  g_bamHeader->text = &l_headerStr[0];
  // if ((g_mnmpOpt->flag & MM_F_OUT_SAM) &&
  //     g_idxReader->n_parts == 1          ) {
  //   // skip header
  // }

  int    l_numSegs = argc - 2;
  char **l_fileName = &argv[2];

  // Construct execution pipeline
  SeqsRead         l_readStg(l_numSegs, l_fileName); 
  MinimapOriginMap l_oriMapStg(l_numThreads);
  MinimapChain     l_chainStg(l_numThreads);
  MinimapAlign     l_alignStg(l_numThreads);
  Reorder          l_reordStg;
  CoordSort        l_sortStg(l_numThreads);
#if 0
  SeqsWrite        l_writeStg(l_numThreads, l_cmdStr.str());
#endif
  SeqsWrite        l_writeStg(l_numThreads);
  MarkDupStage      l_markdupStg(l_numThreads, g_bamHeader);
  BucketSortStage   l_bucketsortStg(g_bamHeader, FLAGS_temp_dir, FLAGS_num_buckets, l_numThreads, FLAGS_compression_level); 


  int l_numStages = 6;
  kestrelFlow::Pipeline l_auxPipe(l_numStages, l_numThreads);
  kestrelFlow::MegaPipe l_mnmpPipe(l_numThreads, 0);

  int l_stg = 0;
  l_auxPipe.addStage(l_stg++, &l_readStg);
  //l_auxPipe.addStage(l_stg++, &l_oriMapStg);
  l_auxPipe.addStage(l_stg++, &l_chainStg);
  l_auxPipe.addStage(l_stg++, &l_alignStg);
  l_auxPipe.addStage(l_stg++, &l_reordStg);
  if (FLAGS_disable_markdup) {
    l_auxPipe.addStage(l_stg++, &l_sortStg);
  }
  else {
    l_auxPipe.addStage(l_stg++, &l_markdupStg);
  }
  if (FLAGS_disable_bucketsort) {
    l_auxPipe.addStage(l_stg++, &l_writeStg);
  }
  else {
    l_auxPipe.addStage(l_stg++, &l_bucketsortStg);
  }
  l_mnmpPipe.addPipeline(&l_auxPipe, 1);


  // Run pipeline
  l_mnmpPipe.start();
  l_mnmpPipe.wait();
  
  l_bucketsortStg.closeFiles();
  
  // Stop timer
  double l_endTime = realtime();

  // Finalize
  std::cerr << "Version: falcon-minimap2 " << VERSION << std::endl;
  std::cerr << "Real time: " << l_endTime - l_startTime << " sec, "
            << "CPUD time: " << cputime() << " sec" << std::endl;

  double sort_start_time = realtime();
  
  if (!FLAGS_disable_bucketsort) {
    kestrelFlow::Pipeline sort_pipeline(4, FLAGS_t);
    
    IndexGenStage     indexgen_stage(
        FLAGS_num_buckets + (!FLAGS_filter_unmap));
    BamReadStage      bamread_stage(FLAGS_temp_dir, g_bamHeader, FLAGS_t);
    BamSortStage      bamsort_stage(FLAGS_t);
    BamWriteStage     bamwrite_stage(
        FLAGS_num_buckets + (!FLAGS_filter_unmap),
        FLAGS_temp_dir, FLAGS_output, g_bamHeader, FLAGS_t);

    sort_pipeline.addStage(0, &indexgen_stage);
    sort_pipeline.addStage(1, &bamread_stage);
    sort_pipeline.addStage(2, &bamsort_stage);
    sort_pipeline.addStage(3, &bamwrite_stage);
    
    kestrelFlow::MegaPipe mp(FLAGS_t, 0); 
    mp.addPipeline(&sort_pipeline, 1); 
  
    mp.start();
    mp.wait();
    
    std::cerr << "sort stage time: " 
      << realtime() - sort_start_time
      << " s" << std::endl;
  }

  // Close index file
  for (mm_idx_t *l_minimizer : g_minimizerNumaList)
    mm_idx_destroy(l_minimizer);
  mm_idx_reader_close(g_idxReader);

  // Clean up
  g_bamHeader->text = NULL;
  g_bamHeader->l_text = 0;
  bam_hdr_destroy(g_bamHeader);
  free(g_mnmpOpt);
  free(g_mnmpIpt);

  return 0;
}
