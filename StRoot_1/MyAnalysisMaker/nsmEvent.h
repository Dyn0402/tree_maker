#ifndef ROOT_nsmEvent
#define ROOT_nsmEvent

#include "TObject.h"

class nsmEvent : public TObject
{
private:
    Float_t vtx_x;
    Float_t vtx_y;
    Float_t vtx_z;
    
    Float_t dca_xy_avg;
    Float_t dca_xy_err;

    UInt_t Nprim;
    UInt_t run;
    UInt_t event_id;

    UInt_t ref2;
    UInt_t ref3;
    
    UInt_t btof;
    
    Float_t qx;
    Float_t qy;

public:
    nsmEvent();
    nsmEvent(Float_t mvtx_x, Float_t mvtx_y, Float_t mvtx_z, Float_t mdca_xy_avg, Float_t mdca_xy_err,
             UInt_t mNprim, UInt_t mrun, UInt_t mevent_id, UInt_t mref2, UInt_t mref3, UInt_t mbtof,
			 Float_t qx, Float_t qy
             );
    
    void SetEventData(Float_t mvtx_x, Float_t mvtx_y, Float_t mvtx_z, Float_t mdca_xy_avg, Float_t mdca_xy_err,
                      UInt_t mNprim, UInt_t mrun, UInt_t mevent_id, UInt_t mref2, UInt_t mref3, UInt_t mbtof,
					  Float_t mqx, Float_t mqy
                      );
    
    Float_t VtxX() { return vtx_x; }
    Float_t VtxY() { return vtx_y; }
    Float_t VtxZ() { return vtx_z; }
    
    Float_t Dca_XY_Avg() { return dca_xy_avg; }
    Float_t Dca_XY_Err() { return dca_xy_err; }

    UInt_t NPrimaries() { return Nprim; }
    UInt_t Run() { return run; }
    UInt_t Event_Id() { return event_id; }
    UInt_t Refmult2() { return ref2; }
    
    UInt_t Btof() { return btof; }

    Float_t Qx() { return qx; }
    Float_t Qy() { return qy; }
    
    ClassDef(nsmEvent,1) //Lambda event class
};

#endif
