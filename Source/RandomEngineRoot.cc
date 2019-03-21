#include "RandomEngineRoot.hh"
#include <iostream>

namespace Garfield {

RandomEngineRoot randomEngine;

RandomEngineRoot::RandomEngineRoot() : RandomEngine(), m_rng(0) {}

RandomEngineRoot::~RandomEngineRoot() {}

void RandomEngineRoot::Seed(const unsigned int s) {
  m_rng.SetSeed(s);
  std::cout << "RandomEngineRoot::Seed:\n"
            << "    Seed: " << m_rng.GetSeed() << "\n";
}

void RandomEngineRoot::Print() {
  std::cout << "RandomEngineRoot::Print:\n"
            << "    Generator type: TRandom3\n"
            << "    Seed: " << m_rng.GetSeed() << "\n";
}

}
