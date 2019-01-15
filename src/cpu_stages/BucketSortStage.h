#ifndef BUCKETSORT
#define BUCKETSORT

#include <cstring>
#include <string.h>
#include <unordered_map>
#include <iostream>
#include <fstream>

#include "kflow/Pipeline.h"
#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"
#include "kflow/SourceStage.h"
#include "kflow/SinkStage.h"
#include "kflow/Common.h"

#include "MnmpData.h"
#include "MnmpCpuStages.h"
#include "MnmpOptions.h"
#include "htslib/sam.h"

#include <boost/atomic.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread.hpp>

//locakable class
class bucketFile :
  public boost::basic_lockable_adapter<boost::mutex>{
  private:
    samFile * fout_;
    bam_hdr_t * head_;
    int32_t id_;
    char* file_path_;
    char* mode_;
    void writeFileHeader();
  public:
    int32_t get_id() {
      return id_;
    }

    bucketFile(bam_hdr_t* head, int32_t id, const char* file_path,
              const char* mode, const htsFormat* fmt): head_(head), id_(id) {
      file_path_ = strdup(file_path);
      mode_ = strdup(mode);
      fout_ = hts_open(file_path_, mode_);
      writeFileHeader();
    }
    ~bucketFile() {
      sam_close(fout_);
      free(file_path_);
      free(mode_);
    }
    void writeFile(std::vector<bam1_t*> vec);
};

class BucketSortStage : public 
    kestrelFlow::MapStage<BamsBatch, int, COMPUTE_DEPTH, 0> {
 public:
  BucketSortStage(bam_hdr_t* head, 
      std::string out_dir, 
      int num_buckets = 1, 
      int n = 1, 
      int l = -1);

  int compute(BamsBatch const & input);

  void closeFiles();

 private:
  int get_bucket_id(bam1_t* read);

  int num_buckets_;
  int bucket_size_;
  bam_hdr_t* head_;
  std::unordered_map<int32_t, bucketFile*> buckets_;
  std::vector<int64_t> accumulate_length_;
  htsFormat fmt_;
};

#endif
