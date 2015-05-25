/*
 * OCP.hpp
 *
 *  Created on: May 16, 2015
 *      Author: zineus
 */

#ifndef OCP_HPP_
#define OCP_HPP_

#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include "ADOL-C_sparseNLP.hpp"
#include "SVector.hpp"
#include "MVector.hpp"

class Phase;
class Guess {
public:
	Guess();
	~Guess();
	MVector nodes;
	MVector x;
	MVector u;
	MVector param;
};

class Config {
public:
	uint max_iter;
	char* hessian_approximation;
	char* NLP_solver;
};

class OCP {
public:
	OCP();
	~OCP();

	Config config;
	Index n_nodes;
	Index n_states;
	Index n_controls;
	Index n_param;
	Index n_events;
	Index n_path;
	Index n_phases;
	Index n_linkages;

	Number *lb_states, *ub_states, *lb_controls, *ub_controls, *lb_param, *ub_param, *lb_path, *ub_path, *lb_events, *ub_events;
	Number lb_t0, ub_t0, lb_tf, ub_tf;

	Guess guess;
	MVector nodes_opt;
	MVector x0_opt;
	MVector xf_opt;
	MVector x_opt;
	MVector u_opt;
	MVector param_opt;
//	char* mgl_marker[5] = 	{"+", "-", "*", "x", "o"};

	void auto_guess_gen();
	ApplicationReturnStatus set_OCP_structure();
	ApplicationReturnStatus NLP_solve();
private:
	void OCPBounds2NLPBounds();
  	SmartPtr<MyADOLC_sparseNLP> myadolc_nlp;
  	SmartPtr<IpoptApplication> app;
  	ApplicationReturnStatus status;
  	SVector sf_u;
  	SVector sf_x;
  	SVector sf_param;
  	SVector t;
};

#endif /* OCP_HPP_ */
