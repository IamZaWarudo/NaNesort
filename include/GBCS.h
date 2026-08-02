#ifndef __GBCS_H__
#define __GBCS_H__

#include <vector>

#include <ddasHit.h>
#include <GPIN.h>
#include <GDSSD.h>
#include <GSSSD.h>
#include <GLaBr.h>
#include <GClover.h>
#include <GPixel.h>


class GBCS {
 public:
  GBCS();
  ~GBCS() = default;

  void Reset();
  void Unpack (const std::vector<ddasHit>& event);

// storing classified events as vectors 
  std::vector<ddasHit> fImplant;
  std::vector<ddasHit> fDecay;
  std::vector<ddasHit> fVeto;

  GPixel fPixel;

// all the detectors 
  GPIN fI2N;  // this only intermitently worked, so largely useless 
  GPIN fI2S;   
  GPIN fI2TAC; // I2S-I2N
  GPIN fPin1;
//  GPIN fPin2;
//  GPIN fPin3;

  GDSSD FrontLowGain;
  GDSSD FrontHighGain;
  GDSSD BackLowGain;
  GDSSD BackHighGain;


  GSSSD fSSSDLow;  // Veto Particles
  GSSSD fSSSDHigh; // gamma spectrum

  GLaBr fLaBr;
  GClover fCrystal;

};




#endif
