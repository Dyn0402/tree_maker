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
class StPicoDstMaker;
class TreeMaker;


void Make_Trees(string input_file_list, string output_dir, int energy) {

	// Load libraries
	cout << "Load" << endl;
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();

	gSystem->Load("StPicoDstMaker");
	gSystem->Load("StPicoEvent");
	gSystem->Load("TreeMaker");

	StChain *chain = new StChain;
	StPicoDstMaker *picoDst_maker = new StPicoDstMaker(2, input_file_list.data());

	// Turn off everything but Primary tracks in order to speed up the analysis and eliminate IO
	picoDst_maker->SetStatus("*", 0);  // Turn off all branches
	picoDst_maker->SetStatus("Event", 1);  // Turn on the Event data (esp. Event number)
	picoDst_maker->SetStatus("Track", 1);  // Turn on the primary track data
	picoDst_maker->SetStatus("BTofPidTraits", 1);
	picoDst_maker->SetStatus("BTofHit", 1);

	picoDst_maker->SetDebug(0);  // Turn off debug information

	TreeMaker *tree_maker = new TreeMaker(picoDst_maker, output_dir, energy);

	int status = chain->Init() ;
	if(status) chain->Fatal(status,"on chain init");

	int num_events = 1e7;
	num_events = picoDst_maker->chain()->GetEntries();

	cout<<"\n############################ Total Event in chain = "<< num_events << " fast entries = " << picoDst_maker->chain()->GetEntriesFast() << "############################\n "<<endl;

	status = 0;
	int event_index = 0;
	int print_period = 1000;
	while(status != 2) {
		chain->Clear();
		if(event_index % print_period == 0) {
			cout << "About to process Event #" << event_index << endl;
		}
		status = chain->Make(event_index++);
		if(status == 2) { --event_index; cout << "Ending on status: " << status << " and " << event_index << " events read " << endl; }
	}

	chain->Finish();

	delete tree_maker;
	delete picoDst_maker;
	delete chain;

	if(event_index > num_events) { cout << endl << endl << "More events found than expected: " << event_index << "/" << num_events << endl << endl; }
	if(event_index < num_events) { cout << endl << endl << "More events found than expected: " << event_index << "/" << num_events << endl << endl; }

	cout << "donzo" << endl;
}
