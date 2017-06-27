// Copyright 2017 <felipeamp>

#ifndef RMQ_H_
#define RMQ_H_

#include <algorithm>
#include <cmath>
#include <memory>
#include <stack>
#include <unordered_map>
#include <vector>


namespace rmq {

template <typename T>
class Naive {
 public:
  explicit Naive(const std::vector<T> &vec) {
    vec_ = std::vector<T>(vec.begin(), vec.end());
  }

  size_t query(size_t start, size_t end) {
    return std::min_element(vec_.cbegin() + start,
                            vec_.cbegin() + end) - vec_.cbegin();
  }

 private:
  std::vector<T> vec_;
};


template <typename T>
class Exact {
public:
  explicit Exact(const std::vector<T> &vec) {
    vec_ = std::vector<T>(vec.begin(), vec.end());

    for (size_t start = 0; start < vec.size(); ++start) {
      std::unordered_map<size_t, size_t> curr;
      for (size_t end = start + 1; end <= vec.size(); ++end) {
        size_t min_index = (std::min_element(vec_.cbegin() + start,
                                             vec_.cbegin() + end)
                            - vec_.cbegin());
        curr.insert({end, min_index});
      }
      exact_rmq_.insert({start, curr});
    }
  }

  size_t query(size_t start, size_t end) {
    return exact_rmq_[start][end];
  }

private:
  std::vector<T> vec_;
  std::unordered_map<size_t, std::unordered_map<size_t, size_t>> exact_rmq_;
};


template <typename T>
class Sparse {
 public:
  explicit Sparse(const std::vector<T> &vec) {
    vec_ = std::vector<T>(vec.begin(), vec.end());

    block_size_ = std::vector<size_t>(vec_.size() + 1, 1);
    size_t curr_index = 1;
    size_t curr_block_size = 1;
    do {
      if (curr_index + curr_block_size > block_size_.size()) {
          std::fill_n(block_size_.begin() + curr_index,
                      block_size_.size() - curr_index,
                      curr_block_size);
      } else {
        std::fill_n(block_size_.begin() + curr_index,
                    curr_block_size,
                    curr_block_size);
      }
      curr_index += curr_block_size;
      curr_block_size *= 2;
    } while (curr_index < block_size_.size());

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
    size_t curr_block_size = block_size_[end - start];
    size_t min_index_1 = block_min_[start][curr_block_size];
    size_t min_index_2 = block_min_[end - curr_block_size][curr_block_size];
    if (vec_[min_index_1] <= vec_[min_index_2])
      return min_index_1;
    else
      return min_index_2;
  }

 private:
  std::vector<std::unordered_map<size_t, size_t>> block_min_;
  std::vector<size_t> block_size_;
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
      block_min_.push_back(
        std::min_element(block_start_it, block_end_it) - vec_.begin());
      high_blocks[block_num] = vec_[block_min_[block_num]];
    }
    block_min_.push_back(
      std::min_element(vec_.begin() + (num_blocks - 1) * block_size_,
                       vec_.end()) - vec_.begin());
    high_blocks[num_blocks - 1] = vec_[block_min_[num_blocks - 1]];

    rmq_high_blocks_ = std::make_unique<rmq::Sparse<T>>(high_blocks);
  }

  size_t query(size_t start, size_t end) {
    size_t start_block_num = start / block_size_;
    size_t end_block_num = end / block_size_;
    if (end_block_num <= start_block_num + 1) {
      return std::min_element(vec_.begin() + start,
                              vec_.begin() + end) - vec_.begin();
    } else {
      size_t min_index_high = rmq_high_blocks_->query(start_block_num + 1,
                                                      end_block_num);
      size_t min_index = block_min_[min_index_high];
      for (size_t index = start;
           index < (start_block_num + 1) * block_size_;
           ++index) {
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


template <typename T>
class FischerHeun {
public:
  explicit FischerHeun(const std::vector<T> &vec) {
    vec_ = std::vector<T>(vec.begin(), vec.end());

    block_size_ = std::floor(0.25 * std::log2(vec_.size()));
    if (block_size_ < 1) {
      block_size_ = 1;
    }
    size_t num_blocks = std::ceil(
      static_cast<double>(vec_.size()) / static_cast<double>(block_size_));

    std::vector<T> high_blocks;
    high_blocks.reserve(num_blocks);
    block_min_.reserve(num_blocks);
    for (size_t block_num = 0; block_num < num_blocks - 1; ++block_num) {
      auto block_start_it = vec_.begin() + block_num * block_size_;
      auto block_end_it = vec_.begin() + (block_num + 1) * block_size_;

      block_min_.push_back(std::min_element(block_start_it, block_end_it)
                           - vec_.begin());

      high_blocks.push_back(vec_[block_min_[block_num]]);

      unsigned long long curr_cartesian_number = get_cartesian_number_(
        vec_.begin() + block_num * block_size_,
        vec_.begin() + (block_num + 1) * block_size_);
      blocks_cartesian_numbers_.push_back(curr_cartesian_number);
      if (rmq_low_blocks_.find(curr_cartesian_number) == rmq_low_blocks_.end()) {
        std::vector<T> block(block_start_it, block_end_it);
        rmq_low_blocks_[curr_cartesian_number] = std::make_unique<
          rmq::Exact<T>>(block);
      }
    }

    // Last block, might be smaller than the others
    block_min_.push_back(
      std::min_element(vec_.begin() + (num_blocks - 1) * block_size_, vec_.end())
                       - vec_.begin());
    high_blocks.push_back(vec_[block_min_[num_blocks - 1]]);

    unsigned long long curr_cartesian_number = get_cartesian_number_(
      vec_.begin() + (num_blocks - 1) * block_size_,
      vec_.end());
    blocks_cartesian_numbers_.push_back(curr_cartesian_number);
    if (rmq_low_blocks_.find(curr_cartesian_number) == rmq_low_blocks_.end()) {
      std::vector<T> block(vec_.begin() + (num_blocks - 1) * block_size_,
                           vec_.end());
      rmq_low_blocks_[curr_cartesian_number] = std::make_unique<
        rmq::Exact<T>>(block);
    }

    rmq_high_blocks_ = std::make_unique<rmq::Sparse<T>>(high_blocks);
  }

  size_t query(size_t start, size_t end) {
    size_t start_block_num = start / block_size_;
    size_t last_block_num = (end - 1) / block_size_;
    if (start_block_num == last_block_num) {
      // Everything contained in one block.
      return (start
              + rmq_low_blocks_[
                  blocks_cartesian_numbers_[start_block_num]]->query(
                    start - start_block_num * block_size_,
                    end - start_block_num * block_size_));
    } else if (start_block_num + 1 == last_block_num) {
      // Everything contained in two blocks.
      size_t first_block_min_index = (
        start
        + rmq_low_blocks_[blocks_cartesian_numbers_[start_block_num]]->query(
          start - start_block_num * block_size_,
          block_size_));
      size_t last_block_min_index = (
        last_block_num * block_size_
        + rmq_low_blocks_[blocks_cartesian_numbers_[last_block_num]]->query(
          0,
          end - last_block_num * block_size_));
      if (vec_[first_block_min_index] <= vec_[last_block_min_index]) {
        return first_block_min_index;
      } else {
        return last_block_min_index;
      }
    } else {
      // Everything contained in 3+ blocks.
      size_t min_index_high = rmq_high_blocks_->query(start_block_num + 1,
                                                      last_block_num);
      size_t high_blocks_min_index = block_min_[min_index_high];
      T high_blocks_min = vec_[high_blocks_min_index];

      size_t first_block_min_index = (
        start
        + rmq_low_blocks_[blocks_cartesian_numbers_[start_block_num]]->query(
          start - start_block_num * block_size_,
          block_size_));
      T first_block_min = vec_[first_block_min_index];

      size_t last_block_min_index = (
        last_block_num * block_size_
        + rmq_low_blocks_[blocks_cartesian_numbers_[last_block_num]]->query(
          0,
          end - last_block_num * block_size_));
      T last_block_min = vec_[last_block_min_index];

      if (high_blocks_min < first_block_min
          && high_blocks_min < last_block_min) {
        return high_blocks_min_index;
      } else if (first_block_min <= last_block_min) {
        return first_block_min_index;
      } else {
        return last_block_min_index;
      }
    }
  }

private:
  size_t block_size_;
  std::unique_ptr<rmq::Sparse<T>> rmq_high_blocks_;
  std::vector<size_t> block_min_;
  std::vector<unsigned long long> blocks_cartesian_numbers_;
  std::unordered_map<unsigned long long,
                     std::unique_ptr<rmq::Exact<T>>> rmq_low_blocks_;
  std::vector<T> vec_;

  unsigned long long get_cartesian_number_(
      const typename std::vector<T>::iterator start_it,
      const typename std::vector<T>::iterator end_it) {
    std::stack<T> s;
    unsigned long long cartesian_number = 0;
    for (auto it = start_it; it < end_it; ++it) {
      T curr_value = *it;
      while (!s.empty() && s.top() >= curr_value) {
        s.pop();
        cartesian_number <<= 1;
      }
      s.push(curr_value);
      cartesian_number <<= 1;
      cartesian_number |= 1ull;
    }
    cartesian_number <<= s.size();
    return cartesian_number;
  }
};

}  // namespace rmq

#endif  // RMQ_H_
