
#include <GLaBr.h>

GLaBr::GLaBr() { Reset();  }
GLaBr::~GLaBr() {}

void GLaBr::Reset() {
  fCharge      = 0;
  fEcal        = 0;
  fTimestamp   = 0;
}


void GLaBr::Unpack(const ddasHit& hit) {
  
  fTimestamp  = hit.GetTimestamp();
  fCharge = hit.GetCharge();
  fEcal   = hit.GetEcal();

}
