#ifndef GCLOVER_H
#define GCLOVER_H

#include <ddasHit.h>
#include <vector>


class GCloverHit {
  public:
    int    fId     {-1};  // crystal 0-63 in hits, clover 0-15 in addbackHits
    double fCharge {0};
    double fEcal   {0};
    double fTime   {-1};

 ClassDef(GCloverHit,1);
};


class GClover {
  public:
  GClover();
  ~GClover();
  
  void Get(int Crystal, const ddasHit& hit);
  void Reset();
  void BuildAddback();

  std::vector<GCloverHit> hits;  
  std::vector<GCloverHit> addbackHits;
  
 ClassDef(GClover,1);
};


#endif
