#include "MnmpCpuStages.h"

#include <gtest/gtest.h>
#include "TestsCommon.h"
#include "TestUtils.h"


TEST_F(CpuStagesTests, ChainAlignStage) {
  // Generate Original Results
  MinimapOriginMap stage_map;
  SeqsBatch i_seqsBatchGolden;
  getSeqsBatch(i_seqsBatchGolden);
  AlignsBatch golden_batch = stage_map.compute(i_seqsBatchGolden);
  
  // Test 
  MinimapChain stage_chain;
  MinimapAlign stage_align;
  SeqsBatch i_seqsBatchTest;
  getSeqsBatch(i_seqsBatchTest);
  ChainsBatch chains_batch = stage_chain.compute(i_seqsBatchTest);
  AlignsBatch test_batch = stage_align.compute(chains_batch);

  ASSERT_EQ(test_batch.m_numSeq, golden_batch.m_numSeq);
  for (int i = 0; i < golden_batch.m_numSeq; i++) {
    ASSERT_EQ(test_batch.m_numReg[i], golden_batch.m_numReg[i]);
    ASSERT_EQ(test_batch.m_repLen[i], golden_batch.m_repLen[i]);
    ASSERT_EQ(test_batch.m_fragGap[i], golden_batch.m_fragGap[i]);
    for (int j = 0; j < golden_batch.m_numReg[i]; j++) {
      ASSERT_EQ(test_batch.m_reg[i][j].cnt, golden_batch.m_reg[i][j].cnt);
      ASSERT_EQ(test_batch.m_reg[i][j].score, golden_batch.m_reg[i][j].score);
    }
  }

  releaseAlignsBatch(golden_batch);
  releaseAlignsBatch(test_batch);
}
