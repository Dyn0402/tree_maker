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
	ParticleVars::ParticleVars();
	ParticleVars::~ParticleVars();

	// Doers
	void clear();

	// Public Members
	vector<float> pt;
	vector<float> phi;
	vector<float> eta;
	vector<float> dca;
	vector<float> nsigma;
	vector<float> beta;
	vector<int> charge;
};


#endif /* PARTICLEVARS_H_ */
