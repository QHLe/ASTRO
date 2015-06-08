/*
 * OCP.hpp
 *
 *  Created on: May 16, 2015
 *      Author: zineus
 */

#ifndef OCP_HPP_
#define OCP_HPP_

enum NLP_SOLVER 	{ma27=0, ma57, ma86, ma97, mumps};
enum OPT_ORDER		{first_order, second_order};

#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include "ADOL-C_sparseNLP.hpp"
#include "SVector.hpp"

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
	Config();
	uint max_iter;
	char* hessian_approximation;
	NLP_SOLVER NLP_solver;
	bool warmstart;
	double NLP_tol;
	OPT_ORDER opt_oder;
};

class OCP {
public:
	OCP();
	~OCP();

	Index n_nodes;
	Index n_states;
	Index n_controls;
	Index n_param;
	Index n_events;
	Index n_path;
	Index n_phases;
	Index n_linkages;

	SVector<double> lb_states, ub_states, lb_controls, ub_controls, lb_param, ub_param, lb_path, ub_path, lb_events, ub_events;
	Number lb_t0, ub_t0, lb_tf, ub_tf;

	Guess 	guess;
	Config 	config;
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
	void set_endpoint_cost(	 double (*)( const  double* ini_states, const  double* fin_states, const  double* param, const  double& t0, const  double& tf, uint phase),
							adouble (*)( const adouble* ini_states, const adouble* fin_states, const adouble* param, const adouble& t0, const adouble& tf, uint phase));
	void set_lagrange_cost(  double (*)( const  double *states, const  double *controls, const  double *param, const  double &time,	uint phase),
							adouble (*)( const adouble *states, const adouble *controls, const adouble *param, const adouble &time,	uint phase));
	void set_derivatives(void (*)( double *states_dot,  double *path, const  double *states, const  double *controls, const  double *param, const  double &time, uint phase),
						 void (*)(adouble *states_dot, adouble *path, const adouble *states, const adouble *controls, const adouble *param, const adouble &time, uint phase));
	void set_events(void (*)( double *events, const  double *ini_states, const  double *fin_states, const  double *param, const  double &t0, const  double &tf, uint phase),
						 void (*)(adouble *events, const adouble *ini_states, const adouble *fin_states, const adouble *param, const adouble &t0, const adouble &tf, uint phase));

private:
	void OCPBounds2NLPBounds();
	void determine_scaling_factors();
  	SmartPtr<MyADOLC_sparseNLP> myadolc_nlp;
  	SmartPtr<IpoptApplication> app;
  	ApplicationReturnStatus status;
  	SVector<double> sf_u;
  	SVector<double> sf_x;
  	SVector<double> sf_path;
  	SVector<double> sf_param;
  	SVector<double> sf_t;
  	SVector<double> sf_events;

//  	SVector sf_f;
//		SVector sf_constraint;

};

#endif /* OCP_HPP_ */
