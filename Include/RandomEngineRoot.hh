#ifndef G_RANDOM_ENGINE_ROOT_H
#define G_RANDOM_ENGINE_ROOT_H

#include <TRandom3.h>

#include "RandomEngine.hh"

namespace Garfield {

/// ROOT random number generator.

class RandomEngineRoot : public RandomEngine {
 public:
  /// Constructor
  RandomEngineRoot();
  /// Destructor
  ~RandomEngineRoot();
  /// Call the random number generator.
  double Draw() override { return m_rng.Rndm(); }
  /// Initialise the random number generator.
  void Seed(const unsigned int s) override;
  /// Print information about the generator used and the seed. 
  void Print() override;

 private:
  TRandom3 m_rng;
};
}

#endif
