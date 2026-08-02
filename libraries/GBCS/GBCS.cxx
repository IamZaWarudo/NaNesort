

#include <ddasHit.h>

#include <GBCS.h>
#include <GPIN.h>
#include <GDSSD.h>
#include <GSSSD.h>
#include <GLaBr.h>
#include <GClover.h>

GBCS::GBCS() { Reset(); }

void GBCS::Reset() {
  fI2N.Reset();
  fI2S.Reset();
  fI2TAC.Reset();
  fPin1.Reset();

  FrontLowGain.Reset();
  FrontHighGain.Reset();
  BackLowGain.Reset();
  BackHighGain.Reset();

  fSSSDLow.Reset();
  fSSSDHigh.Reset();

  fLaBr.Reset();
  fCrystal.Reset();

  fPixel.Reset();
}

// for reference case 0 is strip 1
void GBCS::Unpack(const std::vector<ddasHit>& event) {

  bool Pin1Fired       = false;
  bool frontLowFired   = false;
  bool backLowFired    = false;
  bool frontHighFired  = false;
  bool backHighFired   = false;
  bool VetoIon         = false;

  for(const auto &hit : event) {
  switch(hit.GetId()) {
      case 0 ... 39:     // front High Gain
       FrontHighGain.Get(hit.GetId(), hit);
       frontHighFired  = true; 
        break;
      case 40 ... 79:    // front Low  Gain
       FrontLowGain.Get(hit.GetId() - 40, hit); 
       frontLowFired = true;
        break;
      case 80 ... 119:   // back  High Gain
       BackHighGain.Get(hit.GetId() - 80, hit);
       backHighFired = true; 
        break;
      case 120 ... 159:  // back  Low  Gain
       BackLowGain.Get(hit.GetId() - 120, hit);
       backLowFired  = true; 
        break;
      case 160 ... 175:  // SSSD High Gain strips
       fSSSDHigh.Get(hit.GetId() - 160, hit);
        break;
      case 176:
       fI2N.Get(hit);
        break;
      case 177:
       fI2S.Get(hit); 
        break;
      case 180:
       fI2TAC.Get(hit);
        break;
      case 181:
       fPin1.Get(hit);
       Pin1Fired = true;
        break;
      case 182:
        break;
      case 183:
        break;
      case 192 ... 207:
       fLaBr.Get(hit.GetId() - 192, hit);  // LaBr Crystals
        break;
      case 208 ... 271:  // 16 Clover x 4 crystals
       fCrystal.Get(hit.GetId() - 208, hit);
        break;
      case 272 ... 287:  // SSSD Low Gain strips
       fSSSDLow.Get(hit.GetId() - 272, hit);
       VetoIon = true;
        break;

      default:
      break;
    }
  }
 
 fImplant.clear();
  if (frontLowFired && backLowFired && Pin1Fired && !VetoIon) {
    fImplant = event;}

  fDecay.clear();
  if (frontHighFired && backHighFired && !Pin1Fired && !VetoIon) {
    fDecay = event;}

  fVeto.clear();
  if (VetoIon) {
    fVeto = event;}

  fPixel.BuildImplantPixel(fImplant);   
  fPixel.BuildDecayPixel(fDecay);       
}












