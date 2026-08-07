
#include <GSSSD.h>

GSSSD::GSSSD() { Reset(); }
GSSSD::~GSSSD() {}

void GSSSD::Reset() {
  fId     = 0;
  Strip     = -1;
  fCharge = 0;
  fEcal   = 0;
  fTime   = 0;
}

void GSSSD::Get(int strip, const ddasHit& hit) {

  fId     = hit.GetId();
  fTime   = hit.GetTime();
  fCharge = hit.GetCharge();
  fEcal   = hit.GetEcal();

  int Strip {-1};

  if(fId > 272 && fId < 287){
    Strip = fId - 272;
    }else {
    Strip = fId - 160;
    }

}
