#ifndef __GCORRELATOR_H__
#define __GCORRELATOR_H__

#include <deque>
#include <vector>

#include <GBCS.h>

static constexpr int    fGrid       = 40;
static constexpr double fTicksPerMs = 1.e5;    // 10 nanosecond ticks

struct GCorrPair {
  const GBCS* implant;   // points into fImplant, valid only during the call
  double      dt;        // ms, ALWAYS decay - implant
};


class GCorrelator {
 public:
  GCorrelator();
  ~GCorrelator() = default;

  void AddEvent(const GBCS& bcs);
  void Flush();
  void Print() const;

 private:
  void StoreImplant(const GBCS& imp);
  void StoreDecay(const GBCS& dec);

  void CorrelateDecay(const GBCS& dec);
  int  PruneImplants(double current_time);
  int  FinalizeDecays(double current_time);

  std::deque<GBCS> fImplant;
  std::deque<GBCS> fDecays;

  std::vector<GCorrPair> fCandidates;   // reused every decay

  double fBGWindow    = 1000.0;   // ms, backward
  double fDecayWindow = 3000.0;   // ms, forward
  double fRadius      = 1.0;      // strips, on the centroid

  long fImplantStored  = 0;
  long fImplantsErased  = 0;
  long fDecaysStored    = 0;
  long fDecaysFinalized = 0;

  long fAccepted  = 0;
  long fAmbiguous = 0;
  long fEmpty     = 0;

  long   fOutOfOrder = 0;
  double fTimeCheck  = -1.0;
};


#endif
