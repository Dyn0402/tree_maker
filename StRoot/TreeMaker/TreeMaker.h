/*
 * TreeMaker.h
 *
 *  Created on: Jul 16, 2020
 *      Author: Dylan Neff
 */

#ifndef TREEMAKER_H_
#define TREEMAKER_H_

#include <iostream>
#include <string>

#include "StMaker.h"

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StEvent/StBTofHeader.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoEpdHit.h"

#include "StEvent/StBTofHeader.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TObjArray.h"

#include "BESPars.h"  // pull pars namespace from here
#include "EventVars.h"
#include "ParticleVars.h"

class StMaker;

class StMuDstMaker;
class StMuDst;
class StMuEvent;
class StMuTrack;

class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;

using namespace std;


class TreeMaker : public StMaker {
public:
	// Structors
	TreeMaker(StMuDstMaker* maker);
	TreeMaker(StMuDstMaker* maker, string name, int energy_in, int bes_phase);
	TreeMaker(StPicoDstMaker* maker);
	TreeMaker(StPicoDstMaker* maker, string name, int energy_in, int bes_phase);
	virtual ~TreeMaker();

	// Setters
	void set_energy(int energy_in);  // Set energy being run for cut purposes
	void set_out_file_name(string name);  // Set output root file name

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

	BESPars pars;  // Object with all hardcoded parameters/cuts, need to set_energy_bes(energy, bes_phase)

	string out_file_name;  // Name of output root file
	TFile* out_file;  // Root file to be written
	TTree* tree;  // Tree to be written to out_file

	TH1D* event_cut_hist;
	TH1D* track_cut_hist;
	TH2F* de_dx_pq_hist;
	TH2F* beta_pq_hist;

	// Temp QA plots
//	TH1D* flag_diff_hist;
//	TH1D* nHitsFit_diff_hist;
//	TH1D* nHitsPoss_diff_hist;
//	TH1D* dca_diff_hist;
//	TH2D* dca_prim_glob_hist;
//	TH1D* nHitsFit_diff_post_hist;
//	TH1D* nHitsPoss_diff_post_hist;
//	TH1D* dca_diff_post_hist;
//	TH2D* dca_prim_glob_post_hist;


	int events_read;  // Number of events found and read from input
	int events_processed;  // Number of events processed
	int energy;  // Energy of dataset being read
	int bes_phase;  // Phase of Beam Energy Scan of dataset, 1 (I) or 2 (II)

	EventVars event;
	ParticleVars protons;
	ParticleVars pions;


	ClassDef(TreeMaker, 1)  // Macro for CINT compatibility
};


#endif /* TREEMAKER_H_ */
