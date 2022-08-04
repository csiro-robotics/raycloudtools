// Copyright (c) 2020
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Gavin Catt
#ifndef RAY_RANDOM_H
#define RAY_RANDOM_H

#include <cstdint>
#include <cstdlib>
#include <limits>

#include "raylib/raylibconfig.h"

namespace ray
{
/// Simplistic implementation of the PCG RNG  method http://www.pcg-random.org/.
/// This is about ~5x faster then using the standard std::rand() method with no additional drawbacks.
class RAYLIB_EXPORT PCGRandomGenerator
{
public:
  PCGRandomGenerator()
  {
    seed(12748, 3147, 792751, 14992);
  }  // we can't initialise using std::rand, as it isn't platform independent

  /// The minimum possible result value (inclusive). For compatibility with @c std::shuffle()
  static constexpr unsigned min() { return std::numeric_limits<unsigned int>::min(); }
  /// The maximum possible result value (inclusive). For compatibility with @c std::shuffle()
  static constexpr unsigned max() { return std::numeric_limits<unsigned int>::max(); }

  /// Seed with the given 4 values.
  inline void seed(const unsigned int &a, const unsigned int &b, const unsigned int &c, const unsigned int &d)
  {
    uint64_t s0 = static_cast<uint64_t>(a) << 31 | static_cast<uint64_t>(b);
    uint64_t s1 = static_cast<uint64_t>(c) << 31 | static_cast<uint64_t>(d);

    state_ = 0;
    increment_ = (s1 << 1) | 1;
    (void)operator()();
    state_ += s0;
    (void)operator()();
  }

  /// Return a psudeo random number.
  inline unsigned operator()()
  {
    const uint64_t old_state = state_;
    state_ = old_state * 6364136223846793005ULL + (increment_ | 1);
    const unsigned xor_shifted = static_cast<unsigned>(((old_state >> 18u) ^ old_state) >> 27u);
    const int rot = static_cast<int>(old_state >> 59u);
    return (xor_shifted >> rot) | (xor_shifted << ((-rot) & 31));
  }

  static PCGRandomGenerator &instance();

private:
  /// State of the RNG.
  uint64_t state_ = 0;
  /// Incrementer of the RNG.
  uint64_t increment_ = 0;
};

/// Return a random number using some magic engine behind the hood.
/// Use this over std::rand() as there's a good chance this will be faster/better.
unsigned int RAYLIB_EXPORT rand();
void RAYLIB_EXPORT srand(unsigned int seed);

/// Return a uniformed number between [0,1).
inline double randUniformDouble()
{
  return static_cast<double>(ray::rand() % std::numeric_limits<int>::max()) /
         static_cast<double>(std::numeric_limits<int>::max());
}
}  // namespace ray

#endif  // RAY_RANDOM_H
