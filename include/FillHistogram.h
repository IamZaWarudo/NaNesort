#ifndef __FILLHISTOGRAM_H__
#define __FILLHISTOGRAM_H__

#include <vector>
#include <GBCS.h>
#include <ddasHit.h>


// One-time setup. Call once from main() before the event loop, with the run
// number parsed from the .evt filename. Looks for ../tof/tofNNNN.tof and loads
// it if present; if absent, PID plots fall back to the raw TOF.
void LoadTOF(int run);

// Single source of truth for the time of flight. Every TOF/PID fill goes
// through this so a raw value can never leak into a corrected spectrum.
double CorrectedTOF(const GBCS& bcs);


void FillHistograms(const GBCS& bcs, const std::vector<ddasHit>& hits);

void Channel(const std::vector<ddasHit>& hits);
void PID(const GBCS& bcs);
void Gamma(const GBCS& bcs);
void Other(const GBCS& bcs, const std::vector<ddasHit>& hits);


#endif
