#include <ddasHit.h>

#include <Unpacker.h>
#include <GBCS.h>
#include <GPIN.h>
#include <GDSSD.h>
#include <GSSSD.h>
#include <GLaBr.h>
#include <GClover.h>

Unpacker::Unpacker() { Reset(); }

void Unpacker::Reset() {
  fI2SPin1.Reset();
  fI2SPin2.Reset();
  fI2TAC.Reset();
  fPin1.Reset();
  fPin2.Reset();

  fDSSDLow.Reset();
  fDSSDHigh.Reset();

  fSSSDLow.Reset();
  fSSSDHigh.Reset();

  fLaBr.Reset();
  fCrystal.Reset();
}

// for reference case 0 is strip 1
void Unpacker::Unpack(const std::vector<ddasHit>& event) {

for(const auto &hit : event) {
  switch(hit.GetId()) {
      case 0 ... 39:     // front High Gain
       fDSSDHigh.UnpackFront(hit);
        break;
      case 40 ... 79:    // front Low  Gain
       fDSSDLow.UnpackFront(hit);
        break;
      case 80 ... 119:   // back  High Gain
       fDSSDHigh.UnpackBack(hit);
        break;
      case 120 ... 159:  // back  Low  Gain
       fDSSDLow.UnpackBack(hit);
        break;
      case 160 ... 175:  // SSSD High Gain strips
       fSSSDHigh.Unpack(hit,hit.GetId()-160);
        break;
      case 176:
//     fI2N.Unpack(hit);        // Pin1 - I2N
        break;
      case 177:
       fI2SPin1.Unpack(hit);    // Pin1 - I2S
        break;
      case 179:
//     fI2SPin2.Unpack(hit);    // Pin2 lwk dead 
        break;
      case 180:
       fI2TAC.Unpack(hit);
        break;
      case 181:
       fPin1.Unpack(hit);
        break;
      case 182:
       fPin2.Unpack(hit);
        break;
      case 183:
    // fPin3.Get(hit);    
        break;
      case 192 ... 207:  // SSSD Low Gain strips
       fSSSDLow.Unpack(hit,hit.GetId()-192);
       break;
      case 208 ... 223:
     //  fLaBr.Get(hit.GetId() - 208, hit);  // LaBr Crystals
        break;
      case 224 ... 287:  // 16 Clover x 4 crystals
       fCrystal.Unpack(hit,hit.GetId() - 224);
        break;

      default:
      break;
    }
  }


}


