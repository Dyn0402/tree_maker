/*
 * EventVars.h
 *
 *  Created on: Jul 16, 2020
 *      Author: Dylan Neff
 */

#ifndef EVENTVARS_H_
#define EVENTVARS_H_


class EventVars {
public:
	// Structors
	EventVars();
	~EventVars();

	// Doers
	void clear();

	// Public Members
	int run_num;
	int event_id;
	short refmult;
	short refmult2;
	short refmult3;
	short btof_match;
	short btof_multi;

	float vx, vy, vz;
	float psi_east, psi_east;
	float dca_xy_avg;
	float dca_xy_err;
};



#endif /* EVENTVARS_H_ */
