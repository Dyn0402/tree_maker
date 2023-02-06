#include "Flattener.h"


// Structors
Flattener::Flattener() {
	phi_file_name = "";
	phi_file = NULL;
	ep_file_name = "";
	ep_file = NULL;
	qa_file = NULL;

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
	qa_file = NULL;

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
	qa_file = NULL;

	cent_bins = vector<int>(10);
	iota(cent_bins.begin(), cent_bins.end(), -1);  // Populate cent bins with -1, 0, 1, 2, ..., 8
	phi_types = { "protons", "non-protons" };
	ep_types = { "west", "east" };
}

Flattener::~Flattener() {
	if (qa_file) {
		qa_file->Write();
		qa_file->Close();
	}
}


// Setters
void Flattener::set_qa(string name) {
	qa_file = new TFile(name.data(), "UPDATE");
	for (string phi_type : phi_types) {
		for (int cent_bin : cent_bins) {
			for (int eta_bin = 0; eta_bin < eta_bins; eta_bin++) {
				phi_original_dists[phi_type][cent_bin].push_back({});
				phi_flat_dists[phi_type][cent_bin].push_back({});
			}
		}
	}
}

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
	cout << "Reading phi Fourier coefficients from file..." << endl;
	read_phi_terms();
	cout << "Read phi Fourier coefficients from file" << endl;
}

// Initialize input phi and event plane coefficient files
void Flattener::init_treemaker() {
	phi_file = new TFile(phi_file_name.data(), "READ");
	ep_file = new TFile(ep_file_name.data(), "READ");
	init_phi_terms();
	cout << "Reading phi Fourier coefficients from file..." << endl;
	read_phi_terms();
	cout << "Read phi Fourier coefficients from file." << endl;
	cout << "Reading event plane Fourier coefficients from file..." << endl;
	read_ep_terms();
	cout << "Reading event plane Fourier coefficients from file." << endl;
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
		if (in_string(sin_cos, "cosine")) {
			phi_cos_terms[phi_type][cent_bin][eta_bin][run_key] = (TProfile*)key->ReadObj();
		}
		else if (in_string(sin_cos, "sine")) {
			phi_sin_terms[phi_type][cent_bin][eta_bin][run_key] = (TProfile*)key->ReadObj();
		}
		
		else {
			cout << "Bad sine/cosine read of phi file! " << sin_cos.size() << " " << sin_cos << "  " << file_name << endl;
		}
	}
}

// Read event plane Fourier coefficient TProfiles from file to memory
void Flattener::read_ep_terms() {
	TKey* key;
	TIter key_list(ep_file->GetListOfKeys());
	while ((key = (TKey*)key_list())) {
		string file_name = (string)key->GetName();
		vector<string> file_name_split = split(file_name, '_');
		if (file_name_split.size() != 7) {
			cout << "Bad phi coef object name read, skipping! " << file_name << endl;
			continue;
		}
		string sin_cos = file_name_split[0];
		string ep_type = file_name_split[2];
		int cent_bin = stoi(file_name_split[4]);
		int run_key = stoi(file_name_split[7]);
		if (in_string(sin_cos, "sine")) {
			ep_sin_terms[ep_type][cent_bin][run_key] = (TProfile*)key->ReadObj();
		}
		else if (in_string(sin_cos, "cosine")) {
			ep_cos_terms[ep_type][cent_bin][run_key] = (TProfile*)key->ReadObj();
		}
		else {
			cout << "Bad sine/cosine read of phi file! " << sin_cos.size() << " " << sin_cos << "  " << file_name << endl;
		}
	}
}

// Fill TProfiles with harmonic terms of phi
void Flattener::calc_phi_terms(string particle_type, int cent_bin, float eta, int run, float phi) {
	int eta_bin = get_eta_bin(eta);
	int run_key = get_run_bin_key(run);
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
}

// Fill TProfles with harmonic terms of psi
void Flattener::calc_ep_terms(string ep_type, int cent_bin, int run, float psi) {
	int run_key = get_run_bin_key(run);
	if (ep_sin_terms[ep_type][cent_bin].count(run_key) < 1) {
		ep_file->cd();
		string sin_name = "sine_terms_" + ep_type + "_cent_" + to_string(cent_bin) + "_runkey_" + to_string(run_key);
		ep_sin_terms[ep_type][cent_bin][run_key] = new TProfile(sin_name.data(), "Sine Terms", n_harmonic_high - n_harmonic_low + 1, n_harmonic_low - 0.5, n_harmonic_high + 0.5);
		string cos_name = "cosine_terms_" + ep_type + "_cent_" + to_string(cent_bin) + "_runkey_" + to_string(run_key);
		ep_cos_terms[ep_type][cent_bin][run_key] = new TProfile(cos_name.data(), "Cosine Terms", n_harmonic_high - n_harmonic_low + 1, n_harmonic_low - 0.5, n_harmonic_high + 0.5);
	}
	for (int n = n_harmonic_low; n <= n_harmonic_high; n++) {
		ep_sin_terms[ep_type][cent_bin][run_key]->Fill(n, sin(n * psi));
		ep_cos_terms[ep_type][cent_bin][run_key]->Fill(n, cos(n * psi));
	}
}

// Get shifted phi given original phi based on Fourier coefficients for specific event/track
float Flattener::get_flat_phi(float phi, string particle_type, int cent_bin, float eta, int run) {
	int eta_bin = get_eta_bin(eta);
	int run_key = get_run_bin_key(run);
	if (!phi_sin_terms[particle_type][cent_bin][eta_bin][run_key]) {  // Hopefully just means specific run doesn't exist
		cout << "No run info " << run_key << endl;
		run_key = phi_sin_terms[particle_type][cent_bin][eta_bin].begin()->first;  // Just use any run
		cout << "Using run " << run_key << " instead" << endl;
	}

	TProfile* sin_terms = phi_sin_terms[particle_type][cent_bin][eta_bin][run_key];
	TProfile* cos_terms = phi_cos_terms[particle_type][cent_bin][eta_bin][run_key];

	float dphi = 0.;
	for (int n = n_harmonic_low; n <= n_harmonic_high; n++) {
		dphi += 2. / n * (cos_terms->GetBinContent(n) * sin(n * phi) - sin_terms->GetBinContent(n) * cos(n * phi));
	}

	float shifted_phi = phi + dphi;
	while (shifted_phi >= 2 * M_PI) { shifted_phi -= 2 * M_PI; }
	while (shifted_phi < 0) { shifted_phi += 2 * M_PI; }
	

	if (qa_file) {
		if (phi_original_dists[particle_type][cent_bin][eta_bin].count(run_key) < 1) {
			qa_file->cd();
			string original_name = "original_phi_" + particle_type + "_cent_" + to_string(cent_bin) + "_eta_bin_" + to_string(eta_bin) + "_runkey_" + to_string(run_key);
			phi_original_dists[particle_type][cent_bin][eta_bin][run_key] = new TH1I(original_name.data(), "Original Phi Distribution", 200, 0, 2 * M_PI);
			string flat_name = "flat_phi_" + particle_type + "_cent_" + to_string(cent_bin) + "_eta_bin_" + to_string(eta_bin) + "_runkey_" + to_string(run_key);
			phi_flat_dists[particle_type][cent_bin][eta_bin][run_key] = new TH1I(flat_name.data(), "Flattened Phi Distribution", 200, 0, 2 * M_PI);
		}
		phi_original_dists[particle_type][cent_bin][eta_bin][run_key]->Fill(phi);
		phi_flat_dists[particle_type][cent_bin][eta_bin][run_key]->Fill(shifted_phi);
	}

	return shifted_phi;
}

// Get shifted phi given original phi based on Fourier coefficients for specific event/track
float Flattener::get_flat_ep(float psi, string ep_type, int cent_bin, int run) {
	int run_key = get_run_bin_key(run);
	if (!ep_sin_terms[ep_type][cent_bin][run_key]) {  // Hopefully just means specific run doesn't exist
		cout << "No run info " << run_key << endl;
		run_key = ep_sin_terms[ep_type][cent_bin].begin()->first;  // Just use any run
		cout << "Using run " << run_key << " instead" << endl;
	}

	TProfile* sin_terms = ep_sin_terms[ep_type][cent_bin][run_key];
	TProfile* cos_terms = ep_cos_terms[ep_type][cent_bin][run_key];

	float dpsi = 0.;
	for (int n = n_harmonic_low; n <= n_harmonic_high; n++) {
		dpsi += 2 / n * (cos_terms->GetBinContent(n) * sin(n * psi) - sin_terms->GetBinContent(n) * cos(n * psi));
	}

	float shifted_psi = psi + dpsi;
	while (shifted_psi >= M_PI) { shifted_psi -= M_PI; }
	while (shifted_psi < 0) { shifted_psi += M_PI; }


	if (qa_file) {
		if (psi_original_dists[ep_type][cent_bin].count(run_key) < 1) {
			qa_file->cd();
			string original_name = "original_psi_" + ep_type + "_cent_" + to_string(cent_bin) + "_runkey_" + to_string(run_key);
			psi_original_dists[ep_type][cent_bin][run_key] = new TH1I(original_name.data(), "Original Psi Distribution", 200, 0, M_PI);
			string flat_name = "flat_psi_" + ep_type + "_cent_" + to_string(cent_bin) + "_runkey_" + to_string(run_key);
			psi_flat_dists[ep_type][cent_bin][run_key] = new TH1I(flat_name.data(), "Flattened Phi Distribution", 200, 0, M_PI);
		}
		psi_original_dists[ep_type][cent_bin][run_key]->Fill(psi);
		psi_flat_dists[ep_type][cent_bin][run_key]->Fill(shifted_psi);
	}

	return shifted_psi;
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
	if (eta == eta_max)  return eta_bins - 1;
	float eta_range = eta_max - eta_min;
	return int((eta - eta_min) / eta_range * eta_bins);
}


// Get map key for run_num. For now just truncate last digit of run
int Flattener::get_run_bin_key(int run_num) {
	return int(run_num / run_mod);
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