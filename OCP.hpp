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
	SVector nodes;
	MVector x;
	MVector u;
	SVector param;
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

	double (*lagrange_cost)(const double* states,
							const double* controls,
							const double* param,
							const double& time,
							uint phase);
	double (*endpoint_cost)(const double* ini_states,
							const double* fin_states,
							const double* param,
							const double& t0,
							const double& tf,
							uint phase);
	void (*dynamics)(	double* derivatives,
						double* path_constraint,
						const double* states,
						const double* controls,
						const double* param,
						const double& time,
						uint phase);
	void (*events)(	double* e,
					const double* ini_states,
					const double* fin_states,
					const double* param,
					const double& t0,
					const double& tf,
					uint phase);

	void auto_guess_gen();
	ApplicationReturnStatus NLP_initialization();
	ApplicationReturnStatus NLP_solve();
private:
	void OCPBounds2NLPBounds();
  	SmartPtr<MyADOLC_sparseNLP> myadolc_nlp;
  	SmartPtr<IpoptApplication> app;
  	ApplicationReturnStatus status;
};

#endif /* OCP_HPP_ */
