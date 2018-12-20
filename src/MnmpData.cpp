#include "MnmpData.h"

fragExtSOA *createFragmentExtensionSOA(int i_numFrag) {
  fragExtSOA *o_fragExtSOA = new fragExtSOA;
  o_fragExtSOA->m_segChainsArr = new mm_seg_t*[i_numFrag];

  return o_fragExtSOA;
}

void deleteFragmentExtensionSOA(fragExtSOA *i_fragExtSOA) {
  delete[] i_fragExtSOA->m_segChainsArr;

  delete i_fragExtSOA;
}
