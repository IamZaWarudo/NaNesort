#include <algorithm>
#include <iterator>
#include <cstdio>

#include <GCorrelator.h>
#include <GHistogramer.h>
#include <globals.h>


GCorrelator::GCorrelator() {
  if(fBGWindow > -kDtMin || fDecayWindow > kDtMax)
    printf(RED "\tGCorrelator: windows exceed the corr/dt axis -- "
               "entries will land in overflow and be lost\n" RESET_COLOR);
}


// ---------------------------------------------------------------
// Called once per built event, from ProcessEvent.
// ---------------------------------------------------------------
void GCorrelator::AddEvent(GBCS& bcs, const Unpacker& evt) {

  // MUST come first. EventType() calls GDSSD::HasGoodPosition(), which is
  // the ONLY place fTimestamp / Xpos / Ypos are ever assigned.
  int type = bcs.EventType();
  if(type != 1 && type != 2) return;              // veto / unidentified

  double t = bcs.DSSD.fTimestamp;                 // raw ticks
  int    x = bcs.DSSD.Xpos;
  int    y = bcs.DSSD.Ypos;

  if(x < 0 || x >= kGrid || y < 0 || y >= kGrid) return;

  if(fLastTime > 0 && t < fLastTime) fOutOfOrder++;
  fLastTime = t;

  // fxE = front-strip energy SUM. Not fEcal, which is overwritten by every
  // strip in UnpackFront/UnpackBack and holds whichever was unpacked last.
  if(type == 1) {

    GImplant imp;
    imp.time   = t;
    imp.x      = x;
    imp.y      = y;
    imp.energy = bcs.DSSD.fxE;

    StoreImplant(imp);

  } else {

    GDecay dec;
    dec.time   = t;
    dec.x      = x;
    dec.y      = y;
    dec.energy = bcs.DSSD.fxE;

    // Captured NOW. This decay is not correlated until fBGWindow later, by
    // which point the Unpacker has been Reset() thousands of times over.
    // Reset() sets fEcal = 0, so 0 means "no clover hit this event" -- the
    // same test HasHit() makes, without needing a non-const call.
    dec.gamma = evt.fCrystal.fEcal;

    StoreDecay(dec);
  }

  fDecaysFinalized += FinalizeDecays(t);   // correlate first,
  fImplantsErased  += PruneImplants(t);    // then discard
}


// ---------------------------------------------------------------
// Storage. Neither of these searches anything.
// ---------------------------------------------------------------
void GCorrelator::StoreImplant(const GImplant& imp) {
  fImplantGrid[imp.x][imp.y].push_back(imp);
  fImplantsStored++;
}

void GCorrelator::StoreDecay(const GDecay& dec) {
  fPendingDecays.push_back(dec);
  fDecaysStored++;
}


// ---------------------------------------------------------------
// Release decays whose backward window has closed.
// ---------------------------------------------------------------
int GCorrelator::FinalizeDecays(double current_time) {
  int handled = 0;

  while(!fPendingDecays.empty()) {
    const GDecay& dec = fPendingDecays.front();

    double waited = (current_time - dec.time) / kTicksPerMs;

    // break, NOT continue: the deque is time-ordered, so if the front is
    // not ripe then nothing behind it is either.
    if(waited < fBGWindow) break;

    CorrelateOneDecay(dec);
    fPendingDecays.pop_front();
    handled++;
  }
  return handled;
}


// ---------------------------------------------------------------
// The one search. 3x3 pixels, both signs of dt, one pass.
// ---------------------------------------------------------------
void GCorrelator::CorrelateOneDecay(const GDecay& dec) {

  const int R = fSearchRadius;

  fFwdCandidates.clear();
  fBwdCandidates.clear();

  // ---------- gather ----------
  for(int ix = std::max(0, dec.x - R); ix <= std::min(kGrid-1, dec.x + R); ++ix) {
    for(int iy = std::max(0, dec.y - R); iy <= std::min(kGrid-1, dec.y + R); ++iy) {

      for(const auto& imp : fImplantGrid[ix][iy]) {

        // ALWAYS decay - implant. Never flip this sign anywhere.
        double dt = (dec.time - imp.time) / kTicksPerMs;

        if(dt < -fBGWindow)    continue;      // implant too far AFTER the decay
        if(dt >  fDecayWindow) continue;      // implant too far BEFORE the decay

        GHistogramer::Get().Fill("corr/dt_allpairs", kDtBins, kDtMin, kDtMax, dt);

        if(dt > 0) fFwdCandidates.push_back(dt);
        else       fBwdCandidates.push_back(dt);   // stays negative
      }
    }
  }

  GHistogramer::Get().Fill("corr/ncand_fwd", 20, 0.0, 20.0,
                           (double)fFwdCandidates.size());

  // ---------- forward: signal + background ----------
  if(fFwdCandidates.size() == 1) {
    fFwdAccepted++;
    double dt = fFwdCandidates.front();

    GHistogramer::Get().Fill("corr/dt", kDtBins, kDtMin, kDtMax, dt);

    if(dec.gamma > 0.0) {
      GHistogramer::Get().Fill("corr/gamma", kGammaBins, 0.0, kGammaMax, dec.gamma);
      GHistogramer::Get().Fill("corr/dt_gated", kDtBins, kDtMin, kDtMax, dt);
    }

  } else if(!fFwdCandidates.empty()) {
    fFwdAmbiguous++;                    // counted, never guessed
  } else {
    fFwdEmpty++;
  }

  // ---------- backward: random background, IDENTICAL rule ----------
  if(fBwdCandidates.size() == 1) {
    fBwdAccepted++;
    double dt = fBwdCandidates.front();

    GHistogramer::Get().Fill("corr/dt", kDtBins, kDtMin, kDtMax, dt);

    if(dec.gamma > 0.0) {
      GHistogramer::Get().Fill("corr/gamma_bwd", kGammaBins, 0.0, kGammaMax, dec.gamma);
      GHistogramer::Get().Fill("corr/dt_gated", kDtBins, kDtMin, kDtMax, dt);
    }

  } else if(!fBwdCandidates.empty()) {
    fBwdAmbiguous++;
  } else {
    fBwdEmpty++;
  }
}


// ---------------------------------------------------------------
// Drop implants no decay still to come could ever match.
// ---------------------------------------------------------------
int GCorrelator::PruneImplants(double current_time) {

  // The SUM of both windows. Either alone silently deletes implants that a
  // still-pending decay could have matched. Windows are ms, current_time is
  // ticks -> multiply here (CorrelateOneDecay divides; opposite direction).
  const double tmin_keep =
      current_time - (fBGWindow + fDecayWindow) * kTicksPerMs;

  int erased = 0;

  for(int x = 0; x < kGrid; ++x) {
    for(int y = 0; y < kGrid; ++y) {

      auto& vec = fImplantGrid[x][y];
      if(vec.empty()) continue;

      // Valid only because each pixel's vector is time-ordered: the stream
      // is ordered and we only ever push_back.
      auto first_keep = std::find_if(vec.begin(), vec.end(),
          [tmin_keep](const GImplant& imp) { return imp.time >= tmin_keep; });

      erased += std::distance(vec.begin(), first_keep);
      vec.erase(vec.begin(), first_keep);
    }
  }
  return erased;
}


// ---------------------------------------------------------------
// EOF: correlate whatever is still pending, ripe or not.
// ---------------------------------------------------------------
void GCorrelator::Flush() {
  while(!fPendingDecays.empty()) {
    CorrelateOneDecay(fPendingDecays.front());
    fPendingDecays.pop_front();
    fDecaysFinalized++;
  }
}


void GCorrelator::Print() const {
  printf("\n\tfwd  accepted %ld  ambiguous %ld  empty %ld\n",
         fFwdAccepted, fFwdAmbiguous, fFwdEmpty);
  printf("\tbwd  accepted %ld  ambiguous %ld  empty %ld\n",
         fBwdAccepted, fBwdAmbiguous, fBwdEmpty);
  printf("\timplants stored %ld  erased %ld\n",
         fImplantsStored, fImplantsErased);

  if(fDecaysFinalized != fDecaysStored)
    printf(RED "\t%ld decays never finalized\n" RESET_COLOR,
           fDecaysStored - fDecaysFinalized);

  if(fOutOfOrder)
    printf(RED "\t%ld events out of time order\n" RESET_COLOR, fOutOfOrder);

  printf("\n");
}
