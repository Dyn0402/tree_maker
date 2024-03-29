/*
 * TreeMaker.h
 *
 *  Created on: Jul 16, 2020
 *      Author: Dylan Neff
 */

#ifndef FlattenerPhiEp_H_
#define FlattenerPhiEp_H_

#include <iostream>
#include <string>

#include "StMaker.h"

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StEvent/StBTofHeader.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoEpdHit.h"

#include "StEvent/StBTofHeader.h"

#include "./StRefMultCorr/CentralityMaker.h"
#include "./StRefMultCorr/StRefMultCorr.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TObjArray.h"
#include "TVector2.h"

#include "BESPars.h"
#include "EventVars.h"
#include "ParticleVars.h"
#include "Flattener.h"

class StMaker;

class StMuDstMaker;
class StMuDst;
class StMuEvent;
class StMuTrack;
class StMuPrimaryVertex;

class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;

using namespace std;


class FlattenerPhiEp : public StMaker {
public:
	// Structors
	FlattenerPhiEp(StMuDstMaker* maker);
	FlattenerPhiEp(StMuDstMaker* maker, string name, int energy_in, int bes_phase, string run_type);
	FlattenerPhiEp(StPicoDstMaker* maker);
	FlattenerPhiEp(StPicoDstMaker* maker, string name, int energy_in, int bes_phase, string run_type);
	virtual ~FlattenerPhiEp();

	// Setters
	void set_energy(int energy_in);  // Set energy being run for cut purposes

	// St Doers
	Int_t Init();  // Initialize analysis tools, done once
	Int_t Make();  // Main analysis, performed for each event
	Int_t Finish();  // Finish analysis, close files, clean up

	// Doers
	bool is_bad_event(StMuEvent *mu_event);
	void track_loop(StMuEvent *mu_event);
	bool is_bad_event(StPicoEvent *pico_event);
	void track_loop(StPicoEvent *pico_event);

private:
	StMuDstMaker *muDst_maker;  // MuDstMaker passed in via constructor
	StMuDst *muDst;
	StPicoDstMaker *picoDst_maker;  // PicoDstMaker passed in via constructor
	StPicoDst *picoDst;
	StRefMultCorr *refmultCorrUtil;

	BESPars pars;  // Object with all hardcoded parameters/cuts, need to set_energy_bes(energy, bes_phase)
	EventVars event;
	Flattener flatten;  // Responsible for all flattening and flattening IO

	string run_type;  // "PhiDist" to get Fourier coefficients for phi or "EpDist" for event plane

	int events_read;  // Number of events found and read from input
	int events_processed;  // Number of events processed
	int energy;  // Energy of dataset being read
	int bes_phase;  // Phase of Beam Energy Scan of dataset, 1 (I) or 2 (II)
	int ref_num;  // Reference multiplicity to use for centrality definition

	ClassDef(FlattenerPhiEp, 1)  // Macro for CINT compatibility
};


#endif /* FlattenerPhiEp_H_ */
