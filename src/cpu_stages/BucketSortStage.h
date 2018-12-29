#ifndef BUCKETSORT
#define BUCKETSORT

#include <cstring>
#include <string.h>
#include <unordered_map>

#include "kflow/Pipeline.h"
#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"
#include "kflow/SourceStage.h"
#include "kflow/SinkStage.h"
#include "kflow/Common.h"

#include "MnmpData.h"
#include "MnmpCpuStages.h"
#include "MnmpOptions.h"

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
              const char* mode): head_(head), id_(id) {
      file_path_ = strdup(file_path);
      mode_ = strdup(mode);
      fout_ = sam_open(file_path_, mode_);
      writeFileHeader();
    }
    ~bucketFile() {
      sam_close(fout_);
      free(file_path_);
      free(mode_);
    }
    void writeFile(std::vector<bam1_t*> vec);
};

class BucketSortStage :
  public kestrelFlow::MapStage<BamsBatch, int, COMPUTE_DEPTH, 0> {
  public:
    BucketSortStage(bam_hdr_t* head, std::string out_dir, int num_buckets = 1, int n = 1):
      kestrelFlow::MapStage<BamsBatch, int, COMPUTE_DEPTH, 0>(n), head_(head) {
        accumulate_length_.push_back(0);
        int64_t acc_len = 0;
        for (int i = 0; i < head_->n_targets; i++) {
          acc_len += head_->target_len[i];
          accumulate_length_.push_back(acc_len);
        }
        for (int i = 0; i < num_buckets; i++) {
          //boost::any var = this->getConst("sam_dir");
          //std::string out_dir = boost::any_cast<std::string>(var);
          std::stringstream ss; 
          ss << out_dir << "/part-" << std::setw(6) << std::setfill('0') << i << ".bam";
          const char *modes[] = {"wb", "wb0", "w"};
          bucketFile* tmp_bucket = new bucketFile(head_, i, ss.str().c_str(), modes[FLAGS_output_flag]);
          buckets_[i] = tmp_bucket;
        }
        bucket_size_ = accumulate_length_[head_->n_targets]/num_buckets;
        if (bucket_size_ == 0) {
          throw "bucket_size_ is 0";
        }
        if (accumulate_length_[head_->n_targets]%num_buckets != 0) {
          bucket_size_ += 1;
        }
        #if 0
        std::stringstream interval_file_path;
        interval_file_path << out_dir << "/intervals.list";
        std::ofstream interval_file(interval_file_path.str().c_str());
        int contig_start_pos = 0;
        int contig_id = 0;
        for (int i = 0; i < num_buckets && contig_id < head_->n_targets; i++) {
          int end = contig_start_pos + bucket_size_;
          while (end > head_->target_len[contig_id]) {
            interval_file << contig_id << "\t" << contig_start_pos 
              << "\t" << head_->target_len[contig_id] << "\t" << i << "\n";
            end = end - head_->target_len[contig_id];
            contig_start_pos = 0;
            contig_id += 1;
            if (contig_id >= head_->n_targets) {
              DLOG(INFO) << "unexpected contig id exceeded.";
              break;
            }
          }
          if (contig_id >= head_->n_targets) {
            break;
          }
          interval_file << contig_id << "\t" << contig_start_pos << "\t"
            << end << "\t" << i << "\n";
          contig_start_pos = end;
        }
        interval_file.close();
        #endif
//DLOG(INFO) << "id of bucket 2 " << buckets_[2]->get_id();
      }
    ~BucketSortStage() {
        for (auto it = buckets_.begin(); it != buckets_.end(); ++it) {
          delete it->second;
        }
    }
    int compute(BamsBatch const & input);
  private:
    bam_hdr_t* head_;
    std::unordered_map<int32_t, bucketFile*> buckets_;
    int64_t bucket_size_;
    std::vector<int64_t> accumulate_length_;
    int get_bucket_id(bam1_t* read);
};

#endif
