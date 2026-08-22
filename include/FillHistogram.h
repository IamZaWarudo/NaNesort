#ifndef __FILLHISTOGRAM_H__
#define __FILLHISTOGRAM_H__

#include <vector>
#include <GBCS.h>
#include <ddasHit.h>


void FillHistograms(const GBCS& bcs, const std::vector<ddasHit>& hits);

void Channel(const std::vector<ddasHit>& hits);
void PID(const GBCS& bcs);
void Gamma(const GBCS& bcs);
void Other(const GBCS& bcs, const std::vector<ddasHit>& hits);


#endif
