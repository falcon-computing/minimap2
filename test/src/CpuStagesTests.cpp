#include "MnmpCpuStages.h"

#include <gtest/gtest.h>
#include "kflow/Pipeline.h"
#include "kflow/Queue.h"
#include "TestsCommon.h"
#include "TestUtils.h"

using namespace  kestrelFlow; 

/*
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
*/

TEST_F(CpuStagesTests, KseqTest) {

  KseqsRead        l_kseqsreadStg(2, g_fn);
  KseqsToBseqs     l_kseqs2bseqsStg(2, 1);
  SeqsRead         l_seqsreadStg(2, g_fn);

  kestrelFlow::Pipeline l_auxPipe(2, 1);
  l_auxPipe.addStage(0, &l_kseqsreadStg);
  l_auxPipe.addStage(1, &l_kseqs2bseqsStg);
  
  kestrelFlow::Pipeline l_auxPipe_base(1, 1);
  l_auxPipe_base.addStage(0, &l_seqsreadStg);

  kestrelFlow::MegaPipe mp1(1);
  kestrelFlow::MegaPipe mp2(1);

  mp1.addPipeline(&l_auxPipe);
  mp2.addPipeline(&l_auxPipe_base);

  mp1.start();
  mp1.finalize();
  mp1.wait();

  mp2.start();
  mp2.finalize();
  mp2.wait();

  typedef Queue<SeqsBatch, 64> QIN;
  QIN* iq_test = static_cast<QIN*>(l_auxPipe.getOutputQueue().get());
  QIN* iq_base = static_cast<QIN*>(l_auxPipe_base.getOutputQueue().get());

  SeqsBatch test;
  iq_test->pop(test);

  SeqsBatch base;
  iq_base->pop(base);

  ASSERT_EQ(test.m_numSeq, base.m_numSeq);
  free(test.m_seqs);
  free(base.m_seqs);
}
