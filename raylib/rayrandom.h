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

namespace ray
{
/// Simplistic implementation of the PCG RNG  method http://www.pcg-random.org/.
/// This is about ~5x faster then using the standard std::rand() method with no additional drawbacks.
class RAY_EXPORT PCGRandomGenerator
{
public:
  using result_type = unsigned;

  void PCGRandomGenerator() { seed(std::rand(), std::rand(), std::rand(), std::rand()); }

  /// The minimum possible result value (inclusive). For compatibility with @c std::shuffle()
  static constexpr unsigned min() { return std::numeric_limits<result_type>::min(); }
  /// The maximum possible result value (inclusive). For compatibility with @c std::shuffle()
  static constexpr unsigned max() { return std::numeric_limits<result_type>::max(); }

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
unsigned RAY_EXPORT rand();
/// Return a random number using some magic engine behind the hood which is < then the supplied value.
/// Use this over std::rand() as there's a good chance this will be faster/better.
inline unsigned randWithMax(const unsigned &max)
{
  return rand() % max;
}
/// Return a uniformed number between [0,1).
inline double randUniformDouble()
{
  return static_cast<double>(ras::rand() % std::numeric_limits<int>::max()) /
         static_cast<double>(std::numeric_limits<int>::max());
}
/// Return a uniformed number between the ranges [min,max).
inline double randUniformDouble(const double &min, const double &max)
{
  return min + ((max - min) * randUniformDouble());
}
}  // namespace ray

#endif // RAY_RANDOM_H
