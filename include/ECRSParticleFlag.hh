// 9/28/2026: Shared particle-name -> integer flag lookup.
//
// Replaces the identical 17-way if/else-if string-comparison chains that were
// duplicated in ECRSSteppingAction and ECRSTrackingAction. The lookup table is
// built once (function-local static, thread-safe in C++11) and is read-only
// afterwards, so it is safe to call from multiple worker threads.
//
#ifndef ECRSParticleFlag_hh
#define ECRSParticleFlag_hh 1

#include "globals.hh"
#include <string>
#include <unordered_map>

inline G4int ECRSParticleFlag(const G4String& particleName)
{
  static const std::unordered_map<std::string, G4int> kFlagMap = {
    {"proton", 1},      {"neutron", 2},      {"mu+", 3},   {"mu-", 4},
    {"e+", 5},          {"e-", 6},           {"gamma", 7}, {"pi+", 8},
    {"pi-", 9},         {"C12", 10},         {"C13", 11},  {"He3", 12},
    {"deutron", 13},    {"N14", 14},         {"anti_proton", 15},
    {"anti_neutron", 16}, {"triton", 17}
  };

  auto it = kFlagMap.find(particleName);
  return (it != kFlagMap.end()) ? it->second : 99;
}

#endif
