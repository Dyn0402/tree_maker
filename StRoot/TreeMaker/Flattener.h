#ifndef FLATTENER_H_
#define FLATTENER_H_

#include <iostream>
#include <string>
#include <numeric>
#include <map>

#include "TFile.h"
#include "TProfile.h"
#include "TObject.h"
#include "TKey.h"

using namespace std;

class Flattener {
public:
	// Structors
	Flattener();
	Flattener(string phi_name);
	Flattener(string phi_name, string ep_name);
	virtual ~Flattener();

	// Setters

	// Doers
	void test();
	void draw_phi_term(string particle_type, int cent_bin, int eta_bin, int run_key);
	void init_phi_flattener();
	void init_ep_flattener();
	void init_treemaker();
	void init_phi_terms();
	void read_phi_terms();
	void read_ep_terms();
	void calc_phi_terms(string particle_type, int cent_bin, int eta_bin, int run_key, float phi);
	void calc_ep_terms(string ep_type, int cent_bin, int run_key, float psi);
	void write_phi();
	void write_ep();
	void close_phi_ep();
	float get_flat_phi(float phi, string particle_type, int cent_bin, int eta_bin, int run_key);
	int get_eta_bin(float eta);
	int get_run_bin_key(int run_num);

private:
	string phi_file_name;  // Name of output root file for phi coefficients
	TFile* phi_file;  // Root file to be written for phi coefficients
	string ep_file_name;  // Name of output root file for event plane coefficients
	TFile* ep_file;  // Root file to be written for event plane coefficients

	float eta_min = -1.0;
	float eta_max = 1.0;
	int eta_bins = 20;
	int n_harmonic_low = 1;
	int n_harmonic_high = 12;
	vector<string> phi_types;
	vector<string> ep_types;
	vector<int> cent_bins;

	map<string, map<int, vector<map<int, TProfile*>>>> phi_sin_terms;  // Sine values of particles [particle_type][centrality][eta_bin][run_key]
	map<string, map<int, vector<map<int, TProfile*>>>> phi_cos_terms;  // Cosine values of particles [particle_type][centrality][eta_bin][run_key]

	map<string, map<int, map<int, TProfile*>>> ep_sin_terms;  // Sine values of event planes [ep_type][centrality][run_key]
	map<string, map<int, map<int, TProfile*>>> ep_cos_terms;  // Cosine values of event planes [ep_type][centrality][run_key]
};


// Other functions
// Coppied from Fluctuation_Lib fio_io 2/3/23
bool in_string(string main, string sub);
bool in_string(string main, vector<string> subs, bool all = false);
vector<string> split(string main, char delim = ' ');
vector<string> split(string main, string delim = " ");


#endif /* FLATTENER_H_ */