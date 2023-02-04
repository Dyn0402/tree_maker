#include "Flattener.h"


// Structors
Flattener::Flattener() {
	phi_file_name = "";
	phi_file = NULL;
	ep_file_name = "";
	ep_file = NULL;

	cent_bins = vector<int>(10);
	iota(cent_bins.begin(), cent_bins.end(), -1);  // Populate cent bins with -1, 0, 1, 2, ..., 8
	phi_types = { "protons", "non-protons" };
	ep_types = { "west", "east" };
}

Flattener::Flattener(string phi_name) {
	phi_file_name = phi_name;
	phi_file = NULL;
	ep_file_name = "";
	ep_file = NULL;

	cent_bins = vector<int>(10);
	iota(cent_bins.begin(), cent_bins.end(), -1);  // Populate cent bins with -1, 0, 1, 2, ..., 8
	phi_types = { "protons", "non-protons" };
	ep_types = { "west", "east" };
}

Flattener::Flattener(string phi_name, string ep_name) {
	phi_file_name = phi_name;
	phi_file = NULL;
	ep_file_name = ep_name;
	ep_file = NULL;

	cent_bins = vector<int>(10);
	iota(cent_bins.begin(), cent_bins.end(), -1);  // Populate cent bins with -1, 0, 1, 2, ..., 8
	phi_types = { "protons", "non-protons" };
	ep_types = { "west", "east" };
}

Flattener::~Flattener() {
	// Nothing
}


// Setters


// Getters


// Doers

void Flattener::test() {
	cout << "Flattener working" << endl;
}

void Flattener::draw_phi_term(string particle_type, int cent_bin, int eta_bin, int run_key) {
	phi_sin_terms[particle_type][cent_bin][eta_bin][run_key]->Draw();
}

// Initialize outfile and sin/cos_terms for getting phi coefficients
void Flattener::init_phi_flattener() {
	phi_file = new TFile(phi_file_name.data(), "UPDATE");
	init_phi_terms();
}


// Initialize input phi coefficient file and output event plane coefficient file
void Flattener::init_ep_flattener() {
	phi_file = new TFile(phi_file_name.data(), "READ");
	ep_file = new TFile(ep_file_name.data(), "UPDATE");
	init_phi_terms();
	init_ep_terms();
}

// Initialize input phi and event plane coefficient files
void Flattener::init_treemaker() {
	phi_file = new TFile(phi_file_name.data(), "READ");
	ep_file = new TFile(ep_file_name.data(), "READ");
	init_phi_terms();
	init_ep_terms();
	read_phi_terms();
}

// Initialize phi_sin/cos_terms
void Flattener::init_phi_terms() {
	for (string phi_type : phi_types) {
		for (int cent_bin : cent_bins) {
			for (int eta_bin = 0; eta_bin < eta_bins; eta_bin++) {
				phi_sin_terms[phi_type][cent_bin].push_back({});
				phi_cos_terms[phi_type][cent_bin].push_back({});
			}
		}
	}
}

// Initialize phi_sin/cos_terms
void Flattener::init_ep_terms() {
	for (string ep_type : ep_types) {
		for (string phi_type : phi_types) {
			for (int cent_bin : cent_bins) {
				for (int eta_bin = 0; eta_bin < eta_bins; eta_bin++) {
					ep_sin_terms[ep_type][phi_type][cent_bin].push_back({});
					ep_cos_terms[ep_type][phi_type][cent_bin].push_back({});
				}
			}
		}
	}
}

// Read phi Fourier coefficient TProfiles from file to memory
void Flattener::read_phi_terms() {
	TKey* key;
	TIter key_list(phi_file->GetListOfKeys());
	while ((key = (TKey*)key_list())) {
		string file_name = (string)key->GetName();
		vector<string> file_name_split = split(file_name, '_');
		if (file_name_split.size() != 10) { 
			cout << "Bad phi coef object name read, skipping! " << file_name << endl;
			continue;
		}
		string sin_cos = file_name_split[0];
		string phi_type = file_name_split[2];
		int cent_bin = stoi(file_name_split[4]);
		int eta_bin = stoi(file_name_split[7]);
		int run_key = stoi(file_name_split[9]);
		if (in_string(sin_cos, "sine")) {
			phi_sin_terms[phi_type][cent_bin][eta_bin][run_key] = (TProfile*)key->ReadObj();
		}
		else if (in_string(sin_cos, "cosine")) {
			phi_cos_terms[phi_type][cent_bin][eta_bin][run_key] = (TProfile*)key->ReadObj();
		}
		else {
			cout << "Bad sine/cosine read of phi file! " << sin_cos.size() << " " << sin_cos << "  " << file_name << endl;
		}
	}
}

// Fill TProfiles with harmonic terms of phi
void Flattener::calc_phi_terms(string particle_type, int cent_bin, int eta_bin, int run_key, float phi) {
	if (phi_sin_terms[particle_type][cent_bin][eta_bin].count(run_key) < 1) {
		phi_file->cd();
		string sin_name = "sine_terms_" + particle_type + "_cent_" + to_string(cent_bin) + "_eta_bin_" + to_string(eta_bin) + "_runkey_" + to_string(run_key);
		phi_sin_terms[particle_type][cent_bin][eta_bin][run_key] = new TProfile(sin_name.data(), "Sine Terms", n_harmonic_high - n_harmonic_low + 1, n_harmonic_low - 0.5, n_harmonic_high + 0.5);
		string cos_name = "cosine_terms_" + particle_type + "_cent_" + to_string(cent_bin) + "_eta_bin_" + to_string(eta_bin) + "_runkey_" + to_string(run_key);
		phi_cos_terms[particle_type][cent_bin][eta_bin][run_key] = new TProfile(cos_name.data(), "Cosine Terms", n_harmonic_high - n_harmonic_low + 1, n_harmonic_low - 0.5, n_harmonic_high + 0.5);
	}
	for (int n = n_harmonic_low; n <= n_harmonic_high; n++) {
		phi_sin_terms[particle_type][cent_bin][eta_bin][run_key]->Fill(n, sin(n * phi));
		phi_cos_terms[particle_type][cent_bin][eta_bin][run_key]->Fill(n, cos(n * phi));
	}
	return;
}


// Write and close phi output file
void Flattener::write_phi() {
	phi_file->Write();
	phi_file->Close();
}

// Write and close event plane output file, close input phi file
void Flattener::write_ep() {
	ep_file->Write();
	ep_file->Close();
	phi_file->Close();
}

// Close input event plane and phi files
void Flattener::close_phi_ep() {
	ep_file->Close();
	phi_file->Close();
}


// Calculate eta bin for given eta value
int Flattener::get_eta_bin(float eta) {
	return 1;
	if (eta == eta_max)  return eta_bins - 1;
	float eta_range = eta_max - eta_min;
	return int((eta - eta_min) / eta_range * eta_bins);
}


// Get map key for run_num. For now just truncate last digit of run
int Flattener::get_run_bin_key(int run_num) {
	return int(run_num / 10);
}



// Check if sub is in main.
bool in_string(string main, string sub) {
	if (main.find(sub) != string::npos) {
		return(true);
	}
	return(false);
}



// Other functions
// Check if any subs are in main. If all is true, check if all subs are in main.
bool in_string(string main, vector<string> subs, bool all) {

	for (string sub : subs) {
		if (in_string(main, sub)) {
			if (!all) {
				return(true);
			}
		}
		else if (all) {
			return(false);
		}
	}

	if (all) {
		return(true);
	}
	else {
		return(false);
	}
}


// Emulation of Python split function. Split string into vector of strings on delim.
vector<string> split(string main, char delim) {
	if (main.size() == 0) { return {}; }

	vector<string> split_strings{ "" };
	for (char x : main) {
		if (x == delim) {
			if (split_strings.back() != "") {
				split_strings.push_back("");
			}
		}
		else {
			split_strings.back() += x;
		}
	}
	return(split_strings);
}


// Emulation of Python split function. Split string into vector of strings on delim.
vector<string> split(string main, string delim) {
	if (delim.size() <= 0) { return {}; }
	size_t index = main.find(delim);

	vector<string> split_strings{};
	do {
		split_strings.push_back(main.substr(0, index));
		main = main.substr(index + delim.size());
		index = main.find(delim);
	} while (index != string::npos);

	if (main.size() > 0) { split_strings.push_back(main); }

	return(split_strings);
}