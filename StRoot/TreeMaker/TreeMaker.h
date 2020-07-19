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

#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMaker.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TObjArray.h"

#include "EventVars.h"
#include "ParticleVars.h"
#include "BES1_Parameters.h"  // pull pars namespace from here

class StMaker;
class StMuDstMaker;
class StMuEvent;
class StMuTrack;

using namespace std;


class TreeMaker : public StMaker {
public:
	// Structors
	TreeMaker(StMuDstMaker* maker);
	TreeMaker(StMuDstMaker* maker, string name, int energy_in);
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

private:
	StMuDstMaker* mudst_maker;  // MuDstMaker passed in via constructor

	string out_file_name;  // Name of output root file
	TFile* out_file;  // Root file to be written
	TTree* tree;  // Tree to be written to out_file

	TH1D* event_cut_hist;
	TH1D* track_cut_hist;
	TH2F* de_dx_pq_hist;
	TH2F* beta_pq_hist;

	// Temp QA plots
	TH1D* flag_diff_hist;
	TH1D* nHitsFit_diff_hist;
	TH1D* nHitsPoss_diff_hist;
	TH1D* dca_diff_hist;

	int events_read;  // Number of events found and read from input
	int events_processed;  // Number of events processed
	int energy;  // Energy of dataset being read

	EventVars event;
	ParticleVars protons;
	ParticleVars pions;


	ClassDef(TreeMaker, 1)  // Macro for CINT compatibility
};


#endif /* TREEMAKER_H_ */
