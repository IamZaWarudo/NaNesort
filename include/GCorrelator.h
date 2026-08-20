#ifndef __GCORRELATOR_H__
#define __GCORRELATOR_H__

#include <vector>
#include <deque>

#include <GBCS.h>
#include <Unpacker.h>

// ---- geometry / time ----
static constexpr int    kGrid       = 40;      // DSSD strips 0..39 inclusive
static constexpr double kTicksPerMs = 1.e5;    // 10 ns clock ticks

// ---- corr/dt axis (hardcoded; must cover the windows below) ----
static constexpr double kDtMin  = -1000.0;     // ms
static constexpr double kDtMax  =  5000.0;     // ms
static constexpr int    kDtBins =  6000;      // 10 ms per bin

// ---- gamma axis ----
static constexpr double kGammaMax  = 4000.0;   // widen if clover cal is not keV
static constexpr int    kGammaBins = 4000;     // 1 unit per bin


struct GImplant {
  double time;                    // raw ticks
  int    x;
  int    y;
  double energy;
};

struct GDecay {
  double time;                    // raw ticks
  int    x;
  int    y;
  double energy;
  double gamma;                   // clover fEcal, 0 = no hit this event
};


class GCorrelator {
 public:
  GCorrelator();
  ~GCorrelator() = default;

  // bcs is non-const: EventType() -> HasGoodPosition() mutates GDSSD
  void AddEvent(GBCS& bcs, const Unpacker& evt);
  void Flush();
  void Print() const;

 private:
  void StoreImplant(const GImplant& imp);
  void StoreDecay(const GDecay& dec);

  int  FinalizeDecays(double current_time);
  int  PruneImplants(double current_time);
  void CorrelateOneDecay(const GDecay& dec);

  std::vector<GImplant> fImplantGrid[kGrid][kGrid];
  std::deque<GDecay>    fPendingDecays;

  std::vector<double> fFwdCandidates;   // signed, dt > 0
  std::vector<double> fBwdCandidates;   // signed, dt < 0

  double fBGWindow     = 1000.0;        // ms, backward
  double fDecayWindow  = 5000.0;        // ms, forward
  int    fSearchRadius = 0;             // 1 -> 3x3 pixels

  long fImplantsStored  = 0;
  long fImplantsErased  = 0;
  long fDecaysStored    = 0;
  long fDecaysFinalized = 0;

  long fFwdAccepted = 0, fFwdAmbiguous = 0, fFwdEmpty = 0;
  long fBwdAccepted = 0, fBwdAmbiguous = 0, fBwdEmpty = 0;

  long   fOutOfOrder = 0;
  double fLastTime   = -1.0;
};

#endif
