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

int BucketSortStage::get_bucket_id(bam1_t* read) {
  int32_t contig_id = read->core.tid;
  int32_t read_pos = read->core.pos;
//DLOG(INFO) << "read_pos " << contig_id;
//DLOG(INFO) << "contig_id " << contig_id;
  int64_t acc_pos = accumulate_length_[contig_id] + read_pos;
//DLOG(INFO) << "acc_pos " << acc_pos;
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
//  free(input.bam_buffer);
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
    if (!(vec[i]->core.flag & UNMAP_FLAG)) {
      sam_write1(fout_, head_, vec[i]);
    }
  }
  return;
}
