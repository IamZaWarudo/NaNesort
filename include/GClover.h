#ifndef GCLOVER_H
#define GCLOVER_H

#include <ddasHit.h>
#include <vector>
#include <string>


class GClover {
 public:
  GClover();
  ~GClover();

  double XTalk[16][4][4] = {};
  double CrystalEcal[16][4] = {};    // store Ecal in array
  bool   ActiveCrystal[16][4] = {};  // check Active crystals 
  double EcalCorrection[16][4] = {}; // Apply XTalk correction
  double fAddback[16] = {};          // store addback as one clover

  void AddBack();
  bool LoadXTalk(const char* filename);

  std::vector<ddasHit> fCloverHits;

  void Reset();
  void Unpack(const ddasHit& hit);
  bool HasHit() const;
  
  int    fId;
  double fCharge;
  double fEcal;
  double fTimestamp;

 ClassDef(GClover,1);
};


#endif
