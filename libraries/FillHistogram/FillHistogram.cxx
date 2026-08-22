#include <FillHistogram.h>
#include <GHistogramer.h>
#include <GBCS.h>
#include <ddasHit.h>

// Fill Histograms

void FillHistograms(const GBCS& bcs, const std::vector<ddasHit>& hits){
  
  Channel(hits);
  PID(bcs);
  Gamma(bcs);
  Other(bcs, hits);

  printf("Histograms Filled Succsfully  ദ്ദി(ᵔᗜᵔ)\n");
}

void Channel(const std::vector<ddasHit>& hits){


}


void PID(const GBCS& bcs){



}


void Gamma(const GBCS& bcs){



}


void Other(const GBCS& bcs, const std::vector<ddasHit>& hits){



}




