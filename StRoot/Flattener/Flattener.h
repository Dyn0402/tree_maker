#ifndef Flattener_H_
#define Flattener_H_

#include <iostream>
#include <string>

#include "TFile.h"

using namespace std;

class Flattener {
public:
	// Structors
	Flattener();
	virtual ~Flattener();

	// Setters

	// Doers
	void test();
	void calc_phi_terms(string particle_type, int cent_bin, int eta_bin, int run_key, float phi);
	int get_eta_bin(float eta);
	int get_run_bin_key(int run_num);

private:
	string phi_file_name;  // Name of output root file for phi coefficients
	TFile* phi_file;  // Root file to be written for phi coefficients
	string ep_file_name;  // Name of output root file for event plane coefficients
	TFile* ep_file;  // Root file to be written for event plane coefficients
};

#endif /* Flattener_H_ */