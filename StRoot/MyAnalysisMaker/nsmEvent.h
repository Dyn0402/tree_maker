#ifndef ROOT_nsmEvent
#define ROOT_nsmEvent

#include "TObject.h"

class nsmEvent : public TObject
{
private:
    Float_t vtx_x;
    Float_t vtx_y;
    Float_t vtx_z;
    
    UInt_t Nprim;
    UInt_t run;
    UInt_t ref2;
    
    UInt_t btof;
    
public:
    nsmEvent();
    nsmEvent(Float_t mvtx_x, Float_t mvtx_y, Float_t mvtx_z,
             UInt_t mNprim, UInt_t mrun, UInt_t mref2, UInt_t mbtof
             );
    
    void SetEventData(Float_t mvtx_x, Float_t mvtx_y, Float_t mvtx_z,
                      UInt_t mNprim, UInt_t mrun, UInt_t mref2, UInt_t mbtof
                      );
    
    Float_t VtxX() { return vtx_x; }
    Float_t VtxY() { return vtx_y; }
    Float_t VtxZ() { return vtx_z; }
    
    UInt_t NPrimaries() { return Nprim; }
    UInt_t Run() { return run; }
    UInt_t Refmult2() { return ref2; }
    
    UInt_t Btof() {return btof;}
    
    ClassDef(nsmEvent,1) //Lambda event class
};

#endif
