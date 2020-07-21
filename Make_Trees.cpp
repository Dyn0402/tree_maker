/*
 * Tree_Maker.cpp
 *
 *  Created on: Jul 19, 2020
 *      Author: Dylan Neff
 */


#include <iostream>
#include <string>

#include "TROOT.h"
#include "TSystem.h"

using namespace std;

class StChain;
class StMuDstMaker;
class TreeMaker;


void Make_Trees(string input_file_list, string output_dir, int energy) {
	int num_files = 1e4;

	// Load libraries
	gROOT->Macro("loadMuDst.C");
	gSystem->Load("TreeMaker");
	gSystem->Load("St_base");
	gSystem->Load("StChain");
	gSystem->Load("StUtilities");
	gSystem->Load("StIOMaker");
	gSystem->Load("StarClassLibrary");
	gSystem->Load("StEvent");
	gSystem->Load("StBTofUtil");

	StChain *chain = new StChain;
	StMuDstMaker *muDst_maker = new StMuDstMaker(0, 0, "", input_file_list, "MuDst", num_files);

	// Turn off everything but Primary tracks in order to speed up the analysis and eliminate IO
	muDst_maker->SetStatus("*", 0);  // Turn off all branches
	muDst_maker->SetStatus("MuEvent", 1);  // Turn on the Event data (esp. Event number)
	muDst_maker->SetStatus("PrimaryTracks", 1);  // Turn on the primary track data
	muDst_maker->SetStatus("GlobalTracks", 1);
	muDst_maker->SetStatus("BTofHeader", 1);
	muDst_maker->SetStatus("BTofHit", 1);

	muDst_maker->SetDebug(0);  // Turn off debug information

	TreeMaker *tree_maker = new TreeMaker(muDst_maker, output_dir, energy);

	int num_events = 1e7;
	num_events = muDstMaker->chain()->GetEntries();

	cout<<"\n############################ Total Event in chain = "<< num_events << " fast entries = " << muDstMaker->chain()->GetEntriesFast() << "############################\n "<<endl;


	int status = chain->Init() ;
	if(status) chain->Fatal(status,"on chain init");

	status = 0;
	int event_index = 0;
	int print_period = 1000;
	while(status != 2) {
		chain->Clear();
		if(event_index % print_period == 0) {
			cout << "About to process Event #" << event_index << endl;
		}
		status = chain->Make(event_index++);
		if(status != 2) { cout << "Ending on status: " << status << " and " << event_index << " events read " << endl; }
	}

	chain->Finish();

	delete chain;
}
