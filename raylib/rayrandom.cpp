// Copyright (c) 2019
// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
// ABN 41 687 119 230
//
// Author: Kazys Stepanas
#include "rayrandom.h"

namespace ray
{
PCGRandomGenerator &PCGRandomGenerator::instance()
{
  static PCGRandomGenerator generator;
  return generator;
}

unsigned rand()
{
  /// Static random engine in use for the c-style functions.
  return PCGRandomGenerator::instance()();
}

void srand(unsigned int seed)
{
  // converting a single seed into 4 randomish seeds is an imprecise science...
  PCGRandomGenerator::instance().seed(seed, (seed * 97) % 14746, (13 + seed * 174) % 63721, (seed * 3791 - 27) % 13963);
}
}  // namespace ray
