#ifndef __FILLHISTOGRAM_H__
#define __FILLHISTOGRAM_H__

#include <vector>
#include <GBCS.h>
#include <ddasHit.h>
#include <GCorrelator.h>

void LoadTOF(int run);
double CorrectedTOF(const GBCS& bcs);

void FillHistograms(const GBCS& bcs, const std::vector<ddasHit>& hits);

void Channel(const std::vector<ddasHit>& hits);
void PID(const GBCS& bcs);
void Gamma(const GBCS& bcs);
void Correlation(const GBCS& implant, const GBCS& decay, double dt);
void CorrelationAllPairs(const GBCS& implant, const GBCS& decay, double dt);
void Other(const GBCS& bcs, const std::vector<ddasHit>& hits);


#endif
