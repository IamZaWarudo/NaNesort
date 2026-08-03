#ifndef __GPIXEL_H__
#define __GPIXEL_H__

#include <array>
#include <vector>
#include <ddasHit.h>

class GPixel {
 public:
  GPixel();
  ~GPixel();

  void Reset();

  void BuildImplantPixel (const std::vector<ddasHit>& hits);
  void BuildDecayPixel   (const std::vector<ddasHit>& hits);

  bool ImplantValid() const { return ImplantPixel[0] >= 0.0 && ImplantPixel[1] >= 0.0; }
  bool DecayValid()   const { return DecayPixel[0]   >= 0.0 && DecayPixel[1]   >= 0.0; }

  std::array<double,2> ImplantPixel {-1,-1};
  std::array<double,2> DecayPixel   {-1,-1};

  int ImplantFrontStrip {-1};
  int ImplantBackStrip  {-1};
  int DecayFrontStrip   {-1};
  int DecayBackStrip    {-1};

  double ImplantTime   {-1};
  double DecayTime     {-1};
  double ImplantEnergy {0};
  double DecayEnergy   {0};

  double ImplantFrontE {0};
  double ImplantBackE  {0};
  double DecayFrontE   {0};
  double DecayBackE    {0};

  int ImplantFrontMult {0};
  int ImplantBackMult  {0};
  int DecayFrontMult   {0};
  int DecayBackMult    {0};

 ClassDef(GPixel, 1);
};

#endif
