#include <cstdio>
#include <GCorrelator.h>
#include <FillHistogram.h>
#include <globals.h>


GCorrelator::GCorrelator() {  // self-check 
 if(fBGWindow < -1000 || fDecayWindow > 3000)
    printf(RED "\tGCorrelator: windows exceed the corr/dt axis -- " 
        "entries will land in overflow and be lost\n" RESET_COLOR);
  }


void GCorrelator::AddEvent(const GBCS& bcs){
  
  int type = bcs.EventType();
  double t = bcs.fDSSD.fTimestamp; // current time
  double X = bcs.fDSSD.Xpos;
  double Y = bcs.fDSSD.Ypos;

  if(type != 1 && type != 2) return;   // type check
  if(X < 0 || X >= fGrid || Y < 0 || Y >= fGrid) return;   //position check
  
  if(fTimeCheck > 0 && t < fTimeCheck){ // time ordering check
    printf(RED "\tWOMP WOMP TIME ORDERING ISSUE!! [%ld]\n" RESET_COLOR, fOutOfOrder);
    fflush(stdout);
    fOutOfOrder++;
  }
  fTimeCheck = t;

  // event sorting
  if(type == 1) {
    StoreImplant(bcs);
  }else
  if(type == 2) {
    StoreDecay(bcs);
  }
  
  // counter to keep track
  fDecaysFinalized   += FinalizeDecays(t);  // correlate
  fImplantsErased    += PruneImplants(t);   // prune
}


void GCorrelator::StoreImplant(const GBCS& imp){
  fImplant.push_back(imp);
  fImplantStored++;
}

void GCorrelator::StoreDecay(const GBCS& dec){
  fDecays.push_back(dec);
  fDecaysStored++;
}

int GCorrelator::FinalizeDecays(double current_time) {
  int handled = 0;

  while(!fDecays.empty()){
   const GBCS& dec = fDecays.front();
   double waited = (current_time - dec.fDSSD.fTimestamp) / fTicksPerMs;

   if(waited < fBGWindow) break;

   CorrelateDecay(dec);
   fDecays.pop_front();
   handled++;
  }
  return handled;
}


void GCorrelator::CorrelateDecay(const GBCS& dec) {
  
  const double R = fRadius;
  
  fCandidates.clear();

  for(const auto& imp : fImplant) {
    double dt = (dec.fDSSD.fTimestamp - imp.fDSSD.fTimestamp) / fTicksPerMs;
  
    if(dt < -fBGWindow) break;
    if(dt > fDecayWindow) continue;

    
    double dx = imp.fDSSD.Xpos - dec.fDSSD.Xpos;
    double dy = imp.fDSSD.Ypos - dec.fDSSD.Ypos;
    
    if(dx < -R || dx > R) continue;
    if(dy < -R || dy > R) continue;

    fCandidates.push_back({&imp, dt});

    CorrelationAllPairs(imp, dec, dt);  // using this to see how many pairs am I throwing away
  }

 if(fCandidates.size() != 1){
  if(fCandidates.empty()) fEmpty++;
  else                    fAmbiguous++;
  return;
 }

 fAccepted++;

 const GCorrPair& p = fCandidates.front();
 Correlation(*p.implant, dec, p.dt); // actual correlated pair that passed the check

}

int GCorrelator::PruneImplants(double current_time) {
  const double tmin_keep = current_time - (fBGWindow + fDecayWindow) * fTicksPerMs;

  int erased = 0;
  while(!fImplant.empty() && fImplant.front().fDSSD.fTimestamp < tmin_keep) {
    fImplant.pop_front();
    erased++;
  }
  return erased;
}

void GCorrelator::Flush() {
  while(!fDecays.empty()) {
    CorrelateDecay(fDecays.front());
    fDecays.pop_front();
    fDecaysFinalized++;
  }
}


void GCorrelator::Print() const {
  printf("\n\tcorrelations  accepted %ld  ambiguous %ld  empty %ld\n",
         fAccepted, fAmbiguous, fEmpty);
  printf("\timplants  stored %ld  erased %ld\n",
         fImplantStored, fImplantsErased);
  printf("\tdecays    stored %ld  finalized %ld\n",
         fDecaysStored, fDecaysFinalized);

  if(fOutOfOrder)
    printf(RED "\t%ld events out of time order\n" RESET_COLOR, fOutOfOrder);

  printf("\n");
}

