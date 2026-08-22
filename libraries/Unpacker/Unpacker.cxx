#include <ddasHit.h>

#include <Unpacker.h>
#include <GBCS.h>
#include <GPIN.h>
#include <GDSSD.h>
#include <GSSSD.h>
#include <GLaBr.h>
#include <GClover.h>

Unpacker::Unpacker() {}
Unpacker::~Unpacker() {}



void Unpacker::Unpack(GBCS& bcs, std::vector<ddasHit>& hits) {

bcs.Reset();

for(const auto &hit : hits) {
  switch(hit.GetId()) {
      case 0 ... 39:     // front High Gain
       bcs.fDSSD.UnpackFront(hit);
        break;
      case 40 ... 79:    // front Low  Gain
       bcs.fDSSDLow.UnpackFront(hit);
        break;
      case 80 ... 119:   // back  High Gain
       bcs.fDSSD.UnpackBack(hit);
        break;
      case 120 ... 159:  // back  Low  Gain
       bcs.fDSSDLow.UnpackBack(hit);
        break;
      case 160 ... 175:  // SSSD High Gain strips
       bcs.fSSSDHigh.Unpack(hit,hit.GetId()-160);
        break;
      case 176:
//     fI2N.Unpack(hit);        // Pin1 - I2N
        break;
      case 177:
       bcs.fI2SPin1.Unpack(hit);    // Pin1 - I2S
        break;
      case 179:
//     fI2SPin2.Unpack(hit);    // Pin2 lwk dead 
        break;
      case 180:
       bcs.fI2TAC.Unpack(hit);
        break;
      case 181:
       bcs.fPin1.Unpack(hit);
        break;
      case 182:
       bcs.fPin2.Unpack(hit);
        break;
      case 183:
    // fPin3.Get(hit);    
        break;
      case 192 ... 207:  // SSSD Low Gain strips
       bcs.fSSSDLow.Unpack(hit,hit.GetId()-192);
       break;
      case 208 ... 223:
       //fLaBr.Get(hit.GetId() - 208, hit);  // LaBr Crystals
        break;
      case 224 ... 287:  // 16 Clover x 4 crystals
       bcs.fClover.Unpack(hit);
       bcs.fClover.fCloverHits.push_back(hit);
        break;

      default:
      break;
    }
  }

  bcs.fDSSD.Finalize();
  bcs.fDSSDLow.Finalize();
}


