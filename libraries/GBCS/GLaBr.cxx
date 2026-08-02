
#include <GLaBr.h>

GLaBr::GLaBr() { Reset();  }
GLaBr::~GLaBr() {}

void GLaBr::Reset() {
  fCharge = 0;
  fEcal   = 0;
  fTime   = 0;
}


void GLaBr::Get(int Crystal, const ddasHit& hit) {
  
  fTime   = hit.GetTime();
  fCharge = hit.GetCharge();
  fEcal   = hit.GetEcal();
  // fEcal = hit.GetEcal_labr();

}
