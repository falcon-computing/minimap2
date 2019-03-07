#include "BucketSortStage.h"

#include "MnmpGlobal.h"
#include "MnmpOptions.h"
#include "MnmpWrapper.h"
#include "MnmpUtils.h"
#include "minimap.h"
#include "mmpriv.h"
#include "kalloc.h"

#include "kflow/Common.h"
#include "kflow/Pipeline.h"
#include "kflow/SourceStage.h"
#include "kflow/SinkStage.h"

#define UNMAP_FLAG 4
#define MD_FLAG 1024
#define MAP_IN_PAIR_FLAG 2
#define NOT_PRIM_FLAG 256

std::vector<std::vector<int64_t>> BucketSortStage::get_intervals(int64_t start, int64_t end) {
  if (start < 0 || end < 0) throw("pos less than 0");
  if (start > end) throw("start pos larger than end");
  if (start > accumulate_length_[head_->n_targets]) throw("start pos out of reference boundary");
  std::vector<std::vector<int64_t>> res;
  std::vector<int64_t> tmp;
  // calculate the start contig id;
  int contig_id = 0;
  while (contig_id < head_->n_targets && 
        start >= accumulate_length_[contig_id + 1]) {
    contig_id ++;
  }
  while (contig_id < head_->n_targets &&
        end > accumulate_length_[contig_id + 1]) {
    tmp.push_back(contig_id);
    tmp.push_back(start - accumulate_length_[contig_id]);
    tmp.push_back(accumulate_length_[contig_id + 1] - accumulate_length_[contig_id]);
    res.push_back(tmp);
    tmp.clear();
    start = accumulate_length_[contig_id + 1];
    contig_id ++;
  }
  if (contig_id < head_->n_targets) {
    tmp.push_back(contig_id);
    tmp.push_back(start - accumulate_length_[contig_id]);
    tmp.push_back(end - accumulate_length_[contig_id]);
    res.push_back(tmp);
    tmp.clear();
  }
  return res;
}

int BucketSortStage::bucket_id_calculate(int32_t contig_id, int32_t read_pos) {
  int64_t acc_pos = accumulate_length_[contig_id] + read_pos;
  int large_bucket = (accumulate_length_[head_->n_targets]%num_buckets_)?
                      (accumulate_length_[head_->n_targets]%num_buckets_):
                      num_buckets_;
  int64_t limit = large_bucket * bucket_size_;
  return (acc_pos > limit)?
        (
          (bucket_size_ - 1)?
          (large_bucket + (acc_pos - limit)/(bucket_size_ - 1)):
          (large_bucket)
        ):
        (acc_pos/bucket_size_);
}

int BucketSortStage::get_bucket_id(bam1_t* read) {
  if (read->core.tid == -1) {

    return num_buckets_;
  }
  int32_t contig_id = read->core.tid;
  int32_t read_pos = read->core.pos;
  return bucket_id_calculate(contig_id, read_pos);
}

BucketSortStage::BucketSortStage(
    bam_hdr_t* head, 
    std::string out_dir, 
    int num_buckets, 
    int n, 
    int l): kestrelFlow::MapStage<BamsBatch, int, COMPUTE_DEPTH, 0>(n), 
  num_buckets_(num_buckets), 
  head_(head)
{
  // initialize format
  fmt_.category = sequence_data;
  fmt_.format = bam;
  fmt_.compression = bgzf;
  fmt_.compression_level = l;

  if (!head) {
    throw std::runtime_error("misformat in bam header");
  }

  accumulate_length_.resize(head->n_targets+1);
  accumulate_length_[0] = 0;
  for (int i = 0; i < head_->n_targets; i++) {
    accumulate_length_[i+1] = head_->target_len[i] + accumulate_length_[i];
  }
  int64_t total_length = accumulate_length_[head_->n_targets];
  bucket_size_ = (total_length + num_buckets - 1) / num_buckets;
  int large_bucket = total_length%num_buckets;
  int64_t contig_start_pos = 0;

  const char *modes[] = {"wb", "wb0", "w"};
  // the last bucket is for unmapped reads
  for (int i = 0; i <= num_buckets_; i++) {
    std::stringstream ss; 
    ss << out_dir << "/part-" << std::setw(6) << std::setfill('0') << i;
    std::string bucket_fname = ss.str() + ".bam";
    std::string intv_fname   = ss.str() + ".bed";
    buckets_[i] = new bucketFile(head_, i, 
        bucket_fname.c_str(), 
        modes[FLAGS_output_flag], &fmt_);

    if (i == num_buckets_) break;

    // create interval files
    std::ofstream intv_file(intv_fname.c_str());
    int64_t end = contig_start_pos + bucket_size_ - (int)(i >= large_bucket);
    std::vector<std::vector<int64_t>> intv_vec_vec = get_intervals(contig_start_pos, end);
    for (auto & intv_vec : intv_vec_vec) {
      intv_file << head_->target_name[intv_vec[0]] << "\t"
                << intv_vec[1] << "\t"
                << intv_vec[2] << "\n";
    }
    contig_start_pos = end;
    intv_file.close();
  }
}
  
void BucketSortStage::closeFiles() {
  for (auto it = buckets_.begin(); it != buckets_.end(); ++it) {
    delete it->second;
  }
}

int BucketSortStage::compute(BamsBatch const & input) {
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started BucketSort";
  std::unordered_map<int32_t, std::vector<bam1_t*> > buckets;
  for (int k = 0; k < input.m_numBams; k++) {
        int bucket_id = get_bucket_id(input.m_bams[k]);
        if (buckets.count(bucket_id) != 1) {
          std::vector<bam1_t*> tmp_vec;
          buckets[bucket_id] = tmp_vec;
        }
        buckets[bucket_id].push_back(input.m_bams[k]);
      }
  for (int i = 0; i < buckets_.size(); i++) {
    if (buckets.count(i)) {
      buckets_[i][0].writeFile(buckets[i]);
    }
    for (int j = 0; j < buckets[i].size(); j++) {
      bam_destroy1(buckets[i][j]);
    }
  }
  free(input.m_bams);
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished BucketSort";
  return 0;
}

void bucketFile::writeFileHeader() {
  int status = sam_hdr_write(fout_, head_);
  if (status) {
    ;
  }
  return;
}

void bucketFile::writeFile(std::vector<bam1_t*> vec) {
  boost::lock_guard<bucketFile> guard(*this);
  for (int i = 0; i < vec.size(); i++) {
    if (!(FLAGS_remove_duplicates && (vec[i]->core.flag & MD_FLAG))
        //&& !(vec[i]->core.flag & UNMAP_FLAG)
        //&& (vec[i]->core.flag & MAP_IN_PAIR_FLAG)
        //&& !(vec[i]->core.flag & NOT_PRIM_FLAG)
       ) {
      sam_write1(fout_, head_, vec[i]);
    }
  }
  return;
}
