/*
 * BES1_Parameters.h
 *
 *  Created on: Jul 16, 2020
 *      Author: Dylan Neff
 */

#ifndef BES1_PARAMETERS_H_
#define BES1_PARAMETERS_H_

#include <map>
#include <vector>

using namespace std;


// Hard-coded variables to be defined in BES1_Parameters.cpp

namespace pars {

	extern int branch_buffer;
	extern int branch_split;

	extern map<int, float> vz_max;

	extern map<int, float> vr_max;

	extern map<int, float> vy_offset;

	extern float vertex_min;

	extern map<int, float> vpd_vz_max_diff;

//	extern static const int trig7_arr[];

	extern map<int, vector<int>> triggers;

	extern map<int, vector<int>> bad_runs;
}







#endif /* BES1_PARAMETERS_H_ */
