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

#ifdef BUILD_FPGA
#pragma message "Using FPGA for build"
#ifdef LOCAL_BLAZE
#pragma message "Using local Blaze for build"
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>
#include "blaze/PlatformManager.h"
#include "blaze/AppCommManager.h"
#elif defined(USE_NATIVE)
#pragma message "Using native FPGA routines for build"
#else
#pragma message "Using original C kernels for build"
#endif
#endif

#include "falcon-lic/genome.h"
#include "MnmpGlobal.h"

#include "MnmpUtils.h"
#include "MnmpOptions.h"
#include "MnmpWrapper.h"
#include "MnmpCpuStages.h"
#ifdef BUILD_FPGA
#include "MnmpFpgaStages.h"
#endif

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

  int licret = 0;
  if (0 != (licret = license_verify())) {
    LOG(ERROR) << "Cannot authorize software usage: " << licret;
    LOG(ERROR) << "Please contact support@falcon-computing.com for details.";
    return -1;
  }

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
  kstring_t k_string_tmp = {0,0,0};
  sam_write_rg_line(&k_string_tmp, FLAGS_R.c_str());
  g_bamHeader = sam_hdr_parse(l_headerStr.length(), l_headerStr.c_str());
  g_bamHeader->l_text = l_headerStr.length();
  g_bamHeader->text = &l_headerStr[0];

  int    l_numSegs = argc - 2;
  char **l_fileName = &argv[2];

  // Set up Fpga
#ifdef BUILD_FPGA
#ifdef LOCAL_BLAZE
  // Start Blaze Infrastructure
  blaze::ManagerConf l_blazeConf;
  blaze::PlatformManager *l_blazePlatform = NULL;
  blaze::AppCommManager  *l_blazeComm = NULL;

  do {
    if (!FLAGS_use_fpga)
      break;
    
    int l_confFileHandle = open(FLAGS_blaze_conf.c_str(), O_RDONLY);
    if (l_confFileHandle < 0) {
      LOG(ERROR) << "Cannot find configure file for local blaze: " << FLAGS_blaze_conf;
      LOG(WARNING) << "Not using FPGA acceleration for mmsw";
      FLAGS_use_fpga = false;
      break;
    }
    google::protobuf::io::FileInputStream l_confFile(l_confFileHandle);

    // config manager
    if (!google::protobuf::TextFormat::Parse(&l_confFile, &l_blazeConf)) {
      LOG(ERROR) << "Cannot parse protobuf message";
      LOG(WARNING) << "Not using FPGA acceleration for mmsw";
      FLAGS_use_fpga = false;
      break;
    }
    FLAGS_v = std::max(FLAGS_v, l_blazeConf.verbose());
    // DLOG(INFO) << l_blazeConf.DebugString();
  } while (0);
  
  // start manager
  if (FLAGS_use_fpga) {
    l_blazePlatform = new blaze::PlatformManager(&l_blazeConf);
    l_blazeComm = new blaze::AppCommManager(l_blazePlatform, "127.0.0.1", 1027);

    // Reduce one thread for Blaze manager
    l_numThreads = std::max(l_numThreads-1, 1);
  }
#elif defined(USE_NATIVE)
  // Only check bitstream file in native FPGA mode
  if (FLAGS_use_fpga) {
    boost::filesystem::wpath l_btsm(FLAGS_fpga_path);
    if (FLAGS_use_fpga && !boost::filesystem::exists(l_btsm)) {
      LOG(ERROR) << "Cannot find FPGA bitstream file: " << FLAGS_fpga_path;
      LOG(WARNING) << "Not using FPGA acceleration";
      FLAGS_use_fpga = false;
    }
  }
#endif
  if (!FLAGS_use_fpga && FLAGS_fpga_only) {
    LOG(ERROR) << "Cannot enable CPU or FPGA for computation "
               << "(CPU is disabled due to \"--fpga_only\" option)";
    return -1;
  }
#endif

  l_numThreads = std::max(l_numThreads-FLAGS_extra_threads, 1);

  // Construct execution pipeline
  SeqsRead         l_readStg(l_numSegs, l_fileName); 
  MinimapChain     l_chainStg(l_numThreads);
  MinimapAlign     l_alignStg(l_numThreads);
  Reorder          l_reordStg;
  CoordSort        l_sortStg(l_numThreads);
  SeqsWrite        l_writeStg(l_numThreads);
  MarkDupStage      l_markdupStg(l_numThreads, g_bamHeader); 
  BucketSortStage   l_bucketsortStg(g_bamHeader, FLAGS_temp_dir, FLAGS_num_buckets, l_numThreads, FLAGS_compression_level); 


#ifdef BUILD_FPGA
  MinimapAlignFpga l_alignFpgaStg(FLAGS_fpga_threads);
#endif

  int l_numStages = 6;
  if (!FLAGS_disable_bucketsort) {
    l_numStages += 1;
  }
  kestrelFlow::Pipeline l_auxPipe(l_numStages, l_numThreads);
  kestrelFlow::MegaPipe l_mnmpPipe(l_numThreads, 0);

  int l_stg = 0;
  l_auxPipe.addStage(l_stg++, &l_readStg);
  l_auxPipe.addStage(l_stg++, &l_chainStg);
#ifndef BUILD_FPGA
  l_auxPipe.addStage(l_stg++, &l_alignStg);
#else
  if (!FLAGS_fpga_only)
    l_auxPipe.addStage(l_stg++, &l_alignStg);
  else
    l_auxPipe.addStage(l_stg++, &l_alignFpgaStg);
#endif
  l_auxPipe.addStage(l_stg++, &l_reordStg);
  l_auxPipe.addStage(l_stg++, &l_sortStg);
  if (!FLAGS_disable_markdup) {
    l_auxPipe.addStage(l_stg++, &l_markdupStg);
  }
  if (FLAGS_disable_bucketsort) {
    l_auxPipe.addStage(l_stg++, &l_writeStg);
  }
  else {
    l_auxPipe.addStage(l_stg++, &l_bucketsortStg);
  }
  l_mnmpPipe.addPipeline(&l_auxPipe, 1);

#ifdef BUILD_FPGA
  if (FLAGS_use_fpga && !FLAGS_fpga_only)
    l_auxPipe.addAccxBckStage(2, &l_alignFpgaStg);
#endif

  l_mnmpPipe.addPipeline(&l_auxPipe, 1);

  // Run pipeline
  l_mnmpPipe.start();
  l_mnmpPipe.wait();
  
  l_bucketsortStg.closeFiles();
  
  // Close index file
  for (mm_idx_t *l_minimizer : g_minimizerNumaList)
    mm_idx_destroy(l_minimizer);
  mm_idx_reader_close(g_idxReader);

  // Stop timer
  double l_endTime = realtime();

  // Finalize
  std::cerr << "Version: falcon-minimap2 " << VERSION << std::endl;
  std::cerr << "minimap2 time: " << l_endTime - l_startTime << " sec" << std::endl;

  double sort_start_time = realtime();
  
  if (FLAGS_merge_bams) {
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
    
    std::cerr << "sort time: " 
      << realtime() - sort_start_time
      << " sec" << std::endl;

  }
  // this needs to be performed after BamWriteStage() dealloc
  if (FLAGS_merge_bams) boost::filesystem::remove_all(FLAGS_temp_dir);

  // Clean up
  g_bamHeader->text = NULL;
  g_bamHeader->l_text = 0;
  bam_hdr_destroy(g_bamHeader);
  free(g_mnmpOpt);
  free(g_mnmpIpt);

#ifdef BUILD_FPGA
#ifdef LOCAL_BLAZE
  if (l_blazeComm != NULL)
    delete l_blazeComm;
  if (l_blazePlatform != NULL)
    delete l_blazePlatform; 
#endif
#endif
  return 0;
}
