#ifndef ROOT_nsmTrack
#define ROOT_nsmTrack

#include "TObject.h"

class nsmTrack : public TObject
{
private:

  Float_t pt;
  Float_t p;
  Float_t phi;
  Float_t eta;
  Float_t dca;
  Float_t nsigma;
  Float_t beta;
  Float_t charge;

public:
 nsmTrack();
 nsmTrack(   Float_t mpt, 
	     Float_t mp,
	     Float_t mphi, 
	     Float_t meta, 
	     Float_t mdca,
         Float_t mnsigma,
         Float_t mbeta,
         Float_t mcharge
	     );


 Float_t Pt() { return pt; }
 Float_t P() { return p; }
 Float_t Phi() { return phi; }
 Float_t Eta() { return eta; }
 Float_t Dca() { return dca; }
 Float_t Nsigma() {return nsigma;}
 Float_t Beta() { return beta; }
 Float_t Charge() {return charge;}
 ClassDef(nsmTrack,2) 
};

#endif
