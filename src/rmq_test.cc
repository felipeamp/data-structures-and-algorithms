// Copyright 2017 <felipeamp>

#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

#include "gtest/gtest.h"
#include "rmq.h"

namespace {

constexpr size_t kNumRandomTestsPerSize = 10;
constexpr size_t kNumRandomTestsPerVec = 10;

TEST(RMQ, NaivePreProcessing) {
  std::vector<double> v = {1., 2., 3., 4.};

  rmq::Naive<double> naive(v);
}

TEST(RMQ, NaiveQuery) {
  std::vector<double> v = {1., 2., 3., 4.};

  rmq::Naive<double> naive(v);
  ASSERT_EQ(0, naive.query(0, v.size()));
  ASSERT_EQ(0, naive.query(0, 2));
  ASSERT_EQ(1, naive.query(1, v.size()));
  ASSERT_EQ(2, naive.query(2, v.size()));
}


TEST(RMQ, ExactPreProcessing) {
  std::vector<double> v = {1., 2., 3., 4.};

  rmq::Exact<double> exact(v);
}

TEST(RMQ, ExactQuery) {
  std::vector<double> v = {1., 2., 3., 4.};

  rmq::Exact<double> exact(v);
  ASSERT_EQ(0, exact.query(0, v.size()));
  ASSERT_EQ(0, exact.query(0, 2));
  ASSERT_EQ(1, exact.query(1, v.size()));
  ASSERT_EQ(2, exact.query(2, v.size()));
}

TEST(RMQ, ExactRandomEntries) {
  std::uniform_real_distribution<double> unif_double(0.0, 1.0);
  std::mt19937 re(std::random_device{}());
  auto generator_double = std::bind(unif_double, std::ref(re));
  for (size_t length = 1; length < 100; ++length) {
    std::uniform_int_distribution<size_t> unif_size_t(0, length);
    auto generator_size_t = std::bind(unif_size_t, std::ref(re));

    for (size_t vec_num = 0; vec_num < kNumRandomTestsPerSize; ++vec_num) {
      std::vector<double> vec;
      vec.reserve(length);
      std::generate_n(std::back_inserter(vec),
                      length,
                      std::ref(generator_double));

      for (size_t test_num = 0; test_num < kNumRandomTestsPerVec; ++test_num) {
        size_t start, end;
        do {
          start = generator_size_t();
          end = generator_size_t();
        } while (start == end);
        if (start > end) {
          std::swap(start, end);
        }

        rmq::Naive<double> naive(vec);
        rmq::Exact<double> exact(vec);
        ASSERT_EQ(naive.query(start, end),
                  exact.query(start, end));
      }
    }
  }
}

TEST(RMQ, SparsePreProcessing) {
  std::vector<double> v = {1., 2., 3., 4.};

  rmq::Sparse<double> sparse(v);
}

TEST(RMQ, SparseQuery) {
  std::vector<double> v = {1., 2., 3., 4.};

  rmq::Sparse<double> sparse(v);
  ASSERT_EQ(0, sparse.query(0, v.size()));
  ASSERT_EQ(0, sparse.query(0, 2));
  ASSERT_EQ(1, sparse.query(1, v.size()));
  ASSERT_EQ(2, sparse.query(2, v.size()));
}

TEST(RMQ, SparseRandomEntries) {
  std::uniform_real_distribution<double> unif_double(0.0, 1.0);
  std::mt19937 re(std::random_device{}());
  auto generator_double = std::bind(unif_double, std::ref(re));
  for (size_t length = 1; length < 100; ++length) {
    std::uniform_int_distribution<size_t> unif_size_t(0, length);
    auto generator_size_t = std::bind(unif_size_t, std::ref(re));

    for (size_t vec_num = 0; vec_num < kNumRandomTestsPerSize; ++vec_num) {
      std::vector<double> vec;
      vec.reserve(length);
      std::generate_n(std::back_inserter(vec),
                      length,
                      std::ref(generator_double));

      for (size_t test_num = 0; test_num < kNumRandomTestsPerVec; ++test_num) {
        size_t start, end;
        do {
          start = generator_size_t();
          end = generator_size_t();
        } while (start == end);
        if (start > end) {
          std::swap(start, end);
        }

        rmq::Naive<double> naive(vec);
        rmq::Sparse<double> sparse(vec);
        ASSERT_EQ(naive.query(start, end),
                  sparse.query(start, end));
      }
    }
  }
}

TEST(RMQ, HybridPreProcessing) {
  std::vector<double> v = {1., 2., 3., 4.};

  rmq::Hybrid<double> hybrid(v);
}

TEST(RMQ, HybridQuery) {
  std::vector<double> v = {1., 2., 3., 4.};

  rmq::Hybrid<double> hybrid(v);
  ASSERT_EQ(0, hybrid.query(0, v.size()));
  ASSERT_EQ(0, hybrid.query(0, 2));
  ASSERT_EQ(1, hybrid.query(1, v.size()));
  ASSERT_EQ(2, hybrid.query(2, v.size()));
}

TEST(RMQ, HybridRandomEntries) {
  std::uniform_real_distribution<double> unif_double(0.0, 1.0);
  std::mt19937 re(std::random_device{}());
  auto generator_double = std::bind(unif_double, std::ref(re));
  for (size_t length = 1; length < 100; ++length) {
    std::uniform_int_distribution<size_t> unif_size_t(0, length);
    auto generator_size_t = std::bind(unif_size_t, std::ref(re));

    for (size_t vec_num = 0; vec_num < kNumRandomTestsPerSize; ++vec_num) {
      std::vector<double> vec;
      vec.reserve(length);
      std::generate_n(std::back_inserter(vec),
                      length,
                      std::ref(generator_double));

      for (size_t test_num = 0; test_num < kNumRandomTestsPerVec; ++test_num) {
        size_t start, end;
        do {
          start = generator_size_t();
          end = generator_size_t();
        } while (start == end);
        if (start > end) {
          std::swap(start, end);
        }

        rmq::Naive<double> naive(vec);
        rmq::Hybrid<double> hybrid(vec);
        ASSERT_EQ(naive.query(start, end),
                  hybrid.query(start, end));
      }
    }
  }
}

TEST(RMQ, FischerHeunPreProcessing) {
  std::vector<double> v = {1., 2., 3., 4.};

  rmq::FischerHeun<double> fischer_heun(v);
}

TEST(RMQ, FischerHeunQuery) {
  std::vector<double> v = {1., 2., 3., 4.};

  rmq::FischerHeun<double> fischer_heun(v);
  ASSERT_EQ(0, fischer_heun.query(0, v.size()));
  ASSERT_EQ(0, fischer_heun.query(0, 2));
  ASSERT_EQ(1, fischer_heun.query(1, v.size()));
  ASSERT_EQ(2, fischer_heun.query(2, v.size()));
}

TEST(RMQ, FischerHeunRandomEntries) {
  std::uniform_real_distribution<double> unif_double(0.0, 1.0);
  std::mt19937 re(std::random_device{}());
  auto generator_double = std::bind(unif_double, std::ref(re));
  for (size_t length = 1; length < 100; ++length) {
    std::uniform_int_distribution<size_t> unif_size_t(0, length);
    auto generator_size_t = std::bind(unif_size_t, std::ref(re));

    for (size_t vec_num = 0; vec_num < kNumRandomTestsPerSize; ++vec_num) {
      std::vector<double> vec;
      vec.reserve(length);
      std::generate_n(std::back_inserter(vec),
                      length,
                      std::ref(generator_double));

      for (size_t test_num = 0; test_num < kNumRandomTestsPerVec; ++test_num) {
        size_t start, end;
        do {
          start = generator_size_t();
          end = generator_size_t();
        } while (start == end);
        if (start > end) {
          std::swap(start, end);
        }

        rmq::Naive<double> naive(vec);
        rmq::FischerHeun<double> fischer_heun(vec);
        ASSERT_EQ(naive.query(start, end),
                  fischer_heun.query(start, end));
      }
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
