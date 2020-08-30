/*
 * ParticleVars.h
 *
 *  Created on: Jul 16, 2020
 *      Author: Dylan Neff
 */

#ifndef PARTICLEVARS_H_
#define PARTICLEVARS_H_


#include <vector>

using namespace std;


class ParticleVars {
public:
	// Structors
	ParticleVars();
	~ParticleVars();

	// Doers
	void add_event(float pt, float phi, float eta, float dca, float dca_z, float nsigma, float beta, short charge);
	void clear();

	// Public Members
	vector<float> pt;
	vector<float> phi;
	vector<float> eta;
	vector<float> dca;
	vector<float> dca_z;
	vector<float> nsigma;
	vector<float> beta;
	vector<short> charge;
};


#endif /* PARTICLEVARS_H_ */
