#include "ReorderStage.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>

#include "MnmpGlobal.h"
#include "MnmpOptions.h"
#include "MnmpWrapper.h"
#include "MnmpUtils.h"
#include "minimap.h"
#include "mmpriv.h"
#include "kalloc.h"

#include "htslib/sam.h"


void Reorder::compute(int i_workerId) {
  std::map<int, AlignsBatch> l_inorderCache;
  
  std::vector<AlignsBatch> *l_bundle = new std::vector<AlignsBatch>;
  int l_batchCounter = 0;
  int l_bundleCounter = 0;

  while (true) {
    AlignsBatch i_alignsBatch;
    bool l_ready = this->getInput(i_alignsBatch);
    while (!this->isFinal() && !l_ready) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(10));
      l_ready = this->getInput(i_alignsBatch);
    }
    if (!l_ready) {
      break;
    }

    if (!FLAGS_inorder_output) {
      l_bundle->push_back(i_alignsBatch);
      if (l_bundle->size() >= FLAGS_output_size) {
        AlignsBundle o_alignsBundle;
        o_alignsBundle.m_bundleIdx = l_bundleCounter++;
        o_alignsBundle.m_batches = l_bundle;
        this->pushOutput(o_alignsBundle);
        l_bundle = new std::vector<AlignsBatch>;
      }
    }
    else {
      l_inorderCache[i_alignsBatch.m_batchIdx] = i_alignsBatch;
      while (l_inorderCache.find(l_batchCounter) != l_inorderCache.end()) {
        l_bundle->push_back(l_inorderCache[l_batchCounter]);
        l_inorderCache.erase(l_batchCounter);
        l_batchCounter++;
        if (l_bundle->size() >= FLAGS_output_size) {
          AlignsBundle o_alignsBundle;
          o_alignsBundle.m_bundleIdx = l_bundleCounter++;
          o_alignsBundle.m_batches = l_bundle;
          this->pushOutput(o_alignsBundle);
          l_bundle = new std::vector<AlignsBatch>;
        }
      }
    }
  }

  if (l_bundle->size() > 0) {
    AlignsBundle o_alignsBundle;
    o_alignsBundle.m_bundleIdx = l_bundleCounter++;
    o_alignsBundle.m_batches = l_bundle;
    this->pushOutput(o_alignsBundle);
  }
}
