#include "MnmpGlobal.h"
#include "MnmpOptions.h"
#include "MnmpWrapper.h"
#include "MnmpUtils.h"

#include "BamSortStage.h"
#include "htslib/ksort.h"

static bool bam1_lt(const bam1_t *a, const bam1_t *b) {
  return ((uint64_t)a->core.tid<<32|(a->core.pos+1)<<1|bam_is_rev(a))
       < ((uint64_t)b->core.tid<<32|(b->core.pos+1)<<1|bam_is_rev(b));
}

BamRecord BamSortStage::compute(BamRecord const & input) {
  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Started BamSort()";

  BamRecord output = input;

  if (!FLAGS_disable_sort) {
    //uint64_t start_ts = getUs();
    //sort_bams(output.size, output.bams);
    std::sort(output.bams, output.bams+output.size, bam1_lt);
    //DLOG_IF(INFO, VLOG_IS_ON(1)) << "sorting " << output.size
    //  << " records took " << getUs() - start_ts << " us";
  }

  // write serialized BAM to a buffer first
  // NOTE: not sure how much benefits so far
  //uint64_t start_ts = getUs(); 
  output.fbuf = new BamFileBuffer(0);
  for (int i = 0; i < output.size; i++) {
    if (output.fbuf->write(output.bams[i]) < 0) {
      throw std::runtime_error("cannot convert bam file");
    }
    bam_destroy1(output.bams[i]);
  }
  free(output.bams);

  //DLOG_IF(INFO, VLOG_IS_ON(1)) << "serializing " << output.size
  //  << " records took " << getUs() - start_ts << " us";

  DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finished Bamsort()";
  return output;
}
