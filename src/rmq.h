// Copyright 2017 <felipeamp>

#ifndef RMQ_H_
#define RMQ_H_

#include <algorithm>
#include <cmath>
#include <memory>
#include <unordered_map>
#include <vector>


#include <iostream>

namespace rmq {

template <typename T>
class Naive {
public:
  explicit Naive(const std::vector<T> &vec) {
    vec_ = std::vector<T>(vec.begin(), vec.end());
  }

  size_t query(size_t start, size_t end) {
    return std::min_element(vec_.cbegin() + start, vec_.cbegin() + end) - vec_.cbegin();
  }

private:
  std::vector<T> vec_;
};


template <typename T>
class Sparse {
public:
  explicit Sparse(const std::vector<T> &vec) {
    vec_ = std::vector<T>(vec.begin(), vec.end());

    block_min_ = std::vector<std::unordered_map<size_t, size_t>>(vec_.size());
    block_min_.reserve(vec_.size());
    for (size_t index = 0; index < vec_.size(); ++index) {
      block_min_.emplace_back();
      block_min_[index][1] = index;
    }
    for (size_t block_size = 2; block_size <= vec_.size(); block_size *= 2) {
      for (size_t index = 0; index + block_size <= vec_.size(); ++index) {
        auto min_index_1 = block_min_[index][block_size / 2];
        auto min_index_2 = block_min_[index + block_size / 2][block_size / 2];
        if (vec_[min_index_1] <= vec_[min_index_2]) {
          block_min_[index][block_size] = min_index_1;
        } else {
          block_min_[index][block_size] = min_index_2;
        }
      }
    }
  }

  size_t query(size_t start, size_t end) {
    size_t block_size = (end - start + 1) / 2;
    size_t min_index_1 = block_min_[start][block_size];
    size_t min_index_2 = block_min_[end - block_size][block_size];
    if (vec_[min_index_1] <= vec_[min_index_2])
      return min_index_1;
    else
      return min_index_2;
  }

private:
  std::vector<std::unordered_map<size_t, size_t>> block_min_;
  std::vector<T> vec_;
};

template <typename T>
class Hybrid {
public:
  explicit Hybrid(const std::vector<T> &vec) {
    vec_ = std::vector<T>(vec.begin(), vec.end());

    block_size_ = std::floor(std::log2(vec_.size()));

    size_t num_blocks = std::ceil(
      static_cast<double>(vec_.size()) / static_cast<double>(block_size_));

    std::vector<T> high_blocks(num_blocks);
    block_min_.reserve(num_blocks);
    for (size_t block_num = 0; block_num < num_blocks - 1; ++block_num) {
      auto block_start_it = vec_.begin() + block_num * block_size_;
      auto block_end_it = vec_.begin() + (block_num + 1) * block_size_;
      block_min_.push_back(std::min_element(block_start_it, block_end_it) - vec_.begin());
      high_blocks[block_num] = vec_[block_min_[block_num]];
    }
    block_min_.push_back(
      std::min_element(vec_.begin() + (num_blocks - 1) * block_size_, vec_.end()) - vec_.begin());
    high_blocks[num_blocks - 1] = vec_[block_min_[num_blocks - 1]];

    rmq_high_blocks_ = std::make_unique<rmq::Sparse<T>>(high_blocks);
  }

  size_t query(size_t start, size_t end) {
    size_t start_block_num = start / block_size_;
    size_t end_block_num = end / block_size_;
    if (end_block_num <= start_block_num + 1) {
      return std::min_element(vec_.begin() + start, vec_.begin() + end) - vec_.begin();
    } else {
      size_t min_index_high = rmq_high_blocks_->query(start_block_num + 1, end_block_num);
      size_t min_index = block_min_[min_index_high];
      for (size_t index = start; index < (start_block_num + 1) * block_size_; ++index) {
        if (vec_[index] < vec_[min_index]) {
          min_index = index;
        }
      }
      for (size_t index = end_block_num * block_size_; index < end; ++index) {
        if (vec_[index] < vec_[min_index]) {
          min_index = index;
        }
      }
      return min_index;
    }
  }

private:
  size_t block_size_;
  std::unique_ptr<rmq::Sparse<T>> rmq_high_blocks_;
  std::vector<size_t> block_min_;
  std::vector<T> vec_;
};

// template <typename T>
// class FischerHeun;

}  // namespace rmq

#endif  // RMQ_H_
