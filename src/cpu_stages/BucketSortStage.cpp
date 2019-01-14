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

void BucketSortStage::closeFiles() {
  for (auto it = buckets_.begin(); it != buckets_.end(); ++it) {
    delete it->second;
  }
}

int BucketSortStage::get_bucket_id(bam1_t* read) {
  int32_t contig_id = read->core.tid;
  int32_t read_pos = read->core.pos;
  if (read->core.tid == -1) {
    return num_buckets_;
  }
  int64_t acc_pos = accumulate_length_[contig_id] + read_pos;
  return (acc_pos-1)/bucket_size_;
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
