#ifndef __HEPTAMER_INDEX_HPP__
#define __HEPTAMER_INDEX_HPP__

#include <atomic>
#include <vector>
#include <array>

class HeptamerIndex {
  using AtomicCount = std::atomic<uint64_t>;
public:
  explicit HeptamerIndex();
  std::size_t index(uint64_t heptamer);
private:
  constexpr static uint32_t PossibleHeptamers = 16384;
  const std::array<uint32_t, 7> mult_{{1,4,16,64,256,1024,4096}};
  std::vector<AtomicCount> heptamers_;
};

#endif // __HEPTAMER_INDEX_HPP__
