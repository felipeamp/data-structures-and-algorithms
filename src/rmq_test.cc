// Copyright 2017 <felipeamp>

#include <vector>

#include "gtest/gtest.h"
#include "rmq.h"

namespace {

TEST(RMQ, FirstElem) {
  std::vector<double> v = {1., 2., 3., 4.};

  rmq::Naive<double> naive(v);
  ASSERT_EQ(0, naive.query(0, v.size()));
  ASSERT_EQ(0, naive.query(0, 2));
  ASSERT_EQ(1, naive.query(1, v.size()));
  ASSERT_EQ(2, naive.query(2, v.size()));

  rmq::Sparse<double> sparse(v);
  // ASSERT_EQ(0, sparse.query(0, v.size()));
  // ASSERT_EQ(0, sparse.query(0, 2));
  // ASSERT_EQ(1, sparse.query(1, v.size()));
  // ASSERT_EQ(2, sparse.query(2, v.size()));

  // rmq::Hybrid<double> hybrid(v);
  // ASSERT_EQ(0, hybrid.query(0, v.size()));
  // ASSERT_EQ(0, hybrid.query(0, 2));
  // ASSERT_EQ(1, hybrid.query(1, v.size()));
  // ASSERT_EQ(2, hybrid.query(2, v.size()));

  // rmq::FischerHeun<double> fischer_heun(v);
  // ASSERT_EQ(0, fischer_heun.query(0, v.size()));
  // ASSERT_EQ(0, fischer_heun.query(0, 2));
  // ASSERT_EQ(1, fischer_heun.query(1, v.size()));
  // ASSERT_EQ(2, fischer_heun.query(2, v.size()));
}

}  // namespace

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
