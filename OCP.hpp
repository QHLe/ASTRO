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
	MVector lam_x;
	MVector lam_path;
	MVector lam_events;
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
	bool with_mgl;
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
	Guess 	results;
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
template<class T>
SVector<T>	lin_interpol(const SVector<T> x, const SVector<T> y, const SVector<T> x_new) {
	SVector <T> y_new(x_new.getsize());
	for (uint i = 1; i <= x_new.getsize(); i++) {
		for(uint j = 1; j < x.getsize(); j++) {
/*			if (x_new(i) < x(1) || x(x.getsize()) < x_new(i)) {
				printf("x_new[%d] = %.10e\t	x[%d] = %.10e\n",i,x_new(i),x.getsize(),x(x.getsize()));
				printf("lin_interpot out of range\n");
				exit(1);
			}*/
			if (x(j) <= x_new(i) && x(j+1) >= x_new(i)) {
		//		cout<<i<<"\t"<<j<<endl;
				y_new(i)	= y(j) + (y(j+1) - y(j))/(x(j+1) - x(j))*(x_new(i) - x(j));
				break;
			}
			if (j == 1 || j == x.getsize() - 1) {
				if (x_new(i) < x(1) || x(x.getsize()) < x_new(i)) {
	//				cout<<i<<"\t"<<j<<endl;
					printf("x_new[%d] = %.10e\t	x[%d] = %.10e\n",i,x_new(i),j+1,x(j+1));
					printf("lin_interpot out of range\n");
					//exit(1);
				}
			}

		}
	}
	return y_new;
}

#endif /* OCP_HPP_ */
