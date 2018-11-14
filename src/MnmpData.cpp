#include "MnmpData.h"

fragExtSOA *createFragmentExtensionSOA(int i_numFrag) {
  fragExtSOA *o_fragExtSOA = new fragExtSOA;

  o_fragExtSOA->m_regs0Arr    = new mm_reg1_t*[i_numFrag];
  o_fragExtSOA->m_numRegs0    = new int[i_numFrag];

  o_fragExtSOA->m_hash        = new uint32_t[i_numFrag];

  o_fragExtSOA->m_anchorArr   = new mm128_t*[i_numFrag];

  return o_fragExtSOA;
}

void deleteFragmentExtensionSOA(fragExtSOA *i_fragExtSOA) {
  delete[] i_fragExtSOA->m_regs0Arr;
  delete[] i_fragExtSOA->m_numRegs0;

  delete[] i_fragExtSOA->m_hash;

  delete[] i_fragExtSOA->m_anchorArr;

  delete i_fragExtSOA;
}
