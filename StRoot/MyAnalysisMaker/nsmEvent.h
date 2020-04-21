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
    UInt_t ref3;
    
    UInt_t btof;
    
    Float_t event_plane_ref2;
    Float_t event_plane_ref3;

public:
    nsmEvent();
    nsmEvent(Float_t mvtx_x, Float_t mvtx_y, Float_t mvtx_z,
             UInt_t mNprim, UInt_t mrun, UInt_t mref2, UInt_t mref3, UInt_t mbtof,
			 Float_t event_plane_ref2, Float_t event_plane_ref3
             );
    
    void SetEventData(Float_t mvtx_x, Float_t mvtx_y, Float_t mvtx_z,
                      UInt_t mNprim, UInt_t mrun, UInt_t mref2, UInt mref3, UInt_t mbtof,
					  Float_t mevent_plane_ref2, Float_t event_plane_ref3
                      );
    
    Float_t VtxX() { return vtx_x; }
    Float_t VtxY() { return vtx_y; }
    Float_t VtxZ() { return vtx_z; }
    
    UInt_t NPrimaries() { return Nprim; }
    UInt_t Run() { return run; }
    UInt_t Refmult2() { return ref2; }
    
    UInt_t Btof() { return btof; }

    Float_t EventPlaneRef2() { return event_plane_ref2; }
    Float_t EventPlaneRef3() { return event_plane_ref3; }
    
    ClassDef(nsmEvent,1) //Lambda event class
};

#endif
