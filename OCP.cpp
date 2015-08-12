/*
 * OCP.cpp
 *
 *  Created on: May 16, 2015
 *      Author: zineus
 */

#include <mgl2/mgl.h>
#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include "OCP.hpp"

using namespace Ipopt;

const char* nlp_solver(int enumval) {
	const char* solver_string[] = {"ma27", "ma57", "ma86", "ma97", "mumps"};
	return solver_string[enumval];
}

OCP::OCP() {

	  // Create an instance of your nlp...
	myadolc_nlp = new MyADOLC_sparseNLP();
	  // Create an instance of the IpoptApplication
	app = new IpoptApplication();

	n_nodes			= 0;
	n_states 		= 0;
	n_controls		= 0;
	n_param			= 0;
	n_path			= 0;
	n_events		= 0;
	n_phases		= 0;
	n_linkages		= 0;

	lb_t0 			= 0;
	ub_t0			= 0;

	lb_tf			= 0;
	ub_tf			= 0;

}

OCP::~OCP() {
	// TODO Auto-generated destructor stub
}

ApplicationReturnStatus OCP::set_OCP_structure() {
	Index NLP_n = (n_states + n_controls)*n_nodes + n_param + 2;
	Index NLP_m = ((n_nodes - 1)*n_states + n_events + n_path*n_nodes);

	SVector<uint> structure(8);
	structure(1) = n_phases;
	structure(2) = n_nodes;
	structure(3) = n_states;
	structure(4) = n_controls;
	structure(5) = n_param;
	structure(6) = n_events;
	structure(7) = n_path;
	structure(8) = n_linkages;

	myadolc_nlp->setNLP_structure(NLP_n, NLP_m, structure);

	if (n_states > 0) {
		lb_states.resize(n_states);
		ub_states.resize(n_states);
	}

	if (n_controls > 0) {
		lb_controls.resize(n_controls);
		ub_controls.resize(n_controls);
	}

	if (n_param > 0) {
		lb_param.resize(n_param);
		ub_param.resize(n_param);
	}

	if (n_path > 0) {
		lb_path.resize(n_path);
		ub_path.resize(n_path);
	}

	if (n_events > 0) {
		lb_events.resize(n_events);
		ub_events.resize(n_events);
	}

	app->Options()->SetNumericValue("tol", 1e-6);
	app->Options()->SetStringValue("mu_strategy", "adaptive");
//	app->Options()->SetStringValue("output_file", "ipopt.out");
//	app->Options()->SetStringValue("nlp_scaling_method","gradient-based");
	app->Options()->SetStringValue("linear_solver", "mumps");//	ma86 & ma57 with memmory leakage
	app->Options()->SetStringValue("linear_solver", nlp_solver(config.NLP_solver));//	ma86 & ma57 with memmory leakage
	app->Options()->SetIntegerValue("max_iter", config.max_iter);
	if (config.warmstart) {
		app->Options()->SetStringValue("warm_start_init_point", "yes");
	}

#ifndef SPARSE_HESS
	app->Options()->SetStringValue("hessian_approximation", "limited-memory");
#endif
  	status = app->Initialize();
  	if (status != Solve_Succeeded) {
  		printf("\n\n*** Error during initialization!\n");
  	}
  	return status;
}

ApplicationReturnStatus OCP::NLP_solve() {

	OCPBounds2NLPBounds();

	mglGraph* gr 	= new mglGraph;
	mglData dat_x(n_nodes);
	mglData dat_y(n_nodes);
	mglData x_range(2);
	mglData y_range(2);
	if (config.with_mgl) {
		for (Index i = 0; i < n_nodes; ++i) {
			dat_x[i]	= guess.nodes(i+1,1);
		}

		gr->Box();

		x_range.a[0] = min(guess.nodes);
		x_range.a[1] = max(guess.nodes);

		if (min(guess.x) > 0)
			y_range.a[0] = min(guess.x)*0.9;
		else
			y_range.a[0] = min(guess.x)*1.1;
		if (max(guess.x) < 0)
			y_range.a[1] = max(guess.x)*0.9;
		else
			y_range.a[1] = max(guess.x)*1.1;


		gr->SetRanges(x_range,y_range);
		gr->Axis(); gr->Label('x', "time", 0); gr->Label('y', "x_{guess}", 0);

		for (Index j = 0; j < n_states; ++j) {
			for (Index i = 0; i < n_nodes; ++i) {
				dat_y[i]	= guess.x(i+1,j+1);
			}
			gr->Plot(dat_x,dat_y);
		}

		gr->WriteEPS("guess_x.eps");
		delete gr;

		gr = new mglGraph;
		gr->Box();
		x_range.a[0] = min(guess.nodes);
		x_range.a[1] = max(guess.nodes);

		if (min(guess.u) > 0)
			y_range.a[0] = min(guess.u)*0.9;
		else
			y_range.a[0] = min(guess.u)*1.1;
		if (max(guess.u) < 0)
			y_range.a[1] = max(guess.u)*0.9;
		else
			y_range.a[1] = max(guess.u)*1.1;

		y_range.a[1] = max(guess.u)*1.1;
		gr->SetRanges(x_range,y_range);
		gr->Axis(); gr->Label('x', "time", 0); gr->Label('y', "u_{guess}", 0);

		for (Index j = 0; j < n_controls; ++j) {
			for (Index i = 0; i < n_nodes; ++i) {
				dat_y[i]	= guess.u(i+1,j+1);
			}
			gr->Plot(dat_x,dat_y);
		}


		gr->WriteEPS("guess_u.eps");
		delete gr;
	}

	status = app->OptimizeTNLP(myadolc_nlp);

	if (status == Solve_Succeeded) {
	// Retrieve some statistics about the solve
		Index iter_count = app->Statistics()->IterationCount();
		printf("\n\n*** The problem solved in %d iterations!\n", iter_count);

		Number final_obj = app->Statistics()->FinalObjective();
		printf("\n\n*** The final value of the objective function is %e.\n", final_obj);
	}
	SVector<double> NLP_x 		= myadolc_nlp->get_x_opt();
	SVector<double> NLP_lam 	= myadolc_nlp->get_lam_opt();
	uint NLP_n					= myadolc_nlp->getNLP_n();
//	uint NLP_m					= myadolc_nlp->getNLP_m();
	double t0					= NLP_x(NLP_n - 1);
	double tf 					= NLP_x(NLP_n);
	SVector<double> delta 		= (tf - t0)*myadolc_nlp->getnode_str();

	results.nodes.resize(n_nodes,1);
	results.x.resize(n_nodes, n_states);
	results.u.resize(n_nodes, n_controls);
	results.param.resize(n_param,1);

	results.nodes(1,1) = t0;
	for (Index i = 2; i < n_nodes; i++) {
		results.nodes(i,1) 	= results.nodes(i-1,1) + delta(i);
	}
	results.nodes(n_nodes,1) = tf;
	results.nodes.save("results_nodes.dat");

	Index idx = 1;
	while (idx < n_nodes*(n_states + n_controls) + n_param) {
		for (Index i = 1; i <= n_nodes; i += 1)	{
			for (Index j = 1; j <= n_states; j += 1){
				results.x(i,j) 	= NLP_x(idx);
				idx++;
			}
			for (Index j = 1; j <= n_controls; j += 1){
				results.u(i,j) 	= NLP_x(idx);
				idx++;
			}
		}
		for (Index i = 1; i <= n_param; i += 1)
		{
			results.param(i,1)			= NLP_x(idx);
			idx++;
		}
	}

	results.lam_x.resize(n_nodes - 1, n_states);
	results.lam_path.resize(n_nodes, n_path);
	results.lam_events.resize(n_events, 1);

	Index idx_m = 1;
	for (Index i = 1; i <= n_nodes; i += 1) {
		for (Index j = 1; j <= n_path; j += 1) {
			results.lam_path(i,j)	= -NLP_lam(idx_m);
			idx_m++;
		}
		for (Index j = 1; j <= n_states; j += 1) {
			if(i < n_nodes) {
				results.lam_x(i,j)	= -NLP_lam(idx_m);
				idx_m++;
			}
		}
	}

	for (Index i = 1; i <= n_events; i += 1)	{
		results.lam_events(i,1)		= -NLP_lam(idx_m);
		idx_m++;
	}

	results.x.save("results_x.dat");
	results.u.save("results_u.dat");
	results.param.save("results_param.dat");

	results.lam_x.save("results_lam_x.dat");
	results.lam_path.save("results_lam_path.dat");
	results.lam_events.save("results_lam_events.dat");

	x0_opt.resize(1,n_states);
	xf_opt.resize(1,n_states);
	for (Index i = 1; i <= n_states; i += 1) {
		x0_opt(1,i) 		= results.x(1,i);
		xf_opt(1,i)			= results.x(n_nodes,i);
	}
	x0_opt.save("results_x0.dat");
	xf_opt.save("results_xf.dat");

	if (config.with_mgl) {
		for (Index i = 0; i < n_nodes; ++i) {
			dat_x[i]	= results.nodes(i+1,1);
		}

		gr = new mglGraph;
		gr->Box();

		x_range.a[0] = min(results.nodes);
		x_range.a[1] = max(results.nodes);
		if (min(results.x) > 0)
			y_range.a[0] = min(results.x)*0.9;
		else
			y_range.a[0] = min(results.x)*1.1;
		if (max(results.x) < 0)
			y_range.a[1] = max(results.x)*0.9;
		else
			y_range.a[1] = max(results.x)*1.1;

		gr->SetRanges(x_range,y_range);
		gr->Axis(); gr->Label('x', "time", 0); gr->Label('y', "x_{opt}", 0);

		for (Index j = 0; j < n_states; ++j) {
			for (Index i = 0; i < n_nodes; ++i) {
				dat_y[i]	= results.x(i+1,j+1);
			}
			gr->Plot(dat_x,dat_y);
		}

		gr->WriteEPS("results_x.eps");
		delete gr;

		gr = new mglGraph;
		gr->Box();

		x_range.a[0] = min(results.nodes);
		x_range.a[1] = max(results.nodes);

		cout<<"omg\n";
		if (min(results.u) > 0)
			y_range.a[0] = min(results.u)*0.9;
		else
			y_range.a[0] = min(results.u)*1.1;
		if (max(results.u) < 0)
			y_range.a[1] = max(results.u)*0.9;
		else
			y_range.a[1] = max(results.u)*1.1;

		gr->SetRanges(x_range,y_range);
		gr->Axis(); gr->Label('x', "time", 0); gr->Label('y', "u_{opt}", 0);

		for (Index j = 0; j < n_controls; ++j) {
			for (Index i = 0; i < n_nodes; ++i) {
				dat_y[i]	= results.u(i+1,j+1);
			}
			gr->Plot(dat_x,dat_y);
		}

		gr->WriteEPS("results_u.eps");
	}

	delete gr;
  	return status;
}

void OCP::determine_scaling_factors() {

	sf_x.resize(n_states);
	for (Index i = 1; i <= n_states; i++) {
		if(lb_states(i) != 0 || ub_states(i) != 0) {
			sf_x(i) = max(fabs(lb_states(i)), fabs(ub_states(i)));
		}
		else {
			sf_x(i) = 1;
		}
//		printf("sf_x(%d) = %f\n",i+1,sf_x(i+1));

	}
	sf_u.resize(n_controls);
	for (Index i = 1; i <= n_controls; i++) {
		if(lb_controls(i) != 0 || ub_controls(i) != 0) {
			sf_u(i) = max(fabs(lb_controls(i)), fabs(ub_controls(i)));
		}
		else {
			sf_u(i) = 1;
		}
//		printf("sf_u(%d) = %f\n",i+1,sf_u(i+1));
	}
	sf_path.resize(n_path);
	for (Index i = 1; i <= n_path; i++) {
		if(lb_path(i) != 0 || ub_path(i) != 0) {
			sf_path(i) = max(fabs(lb_path(i)), fabs(ub_path(i)));
		}
		else {
			sf_path(i) = 1;
		}
//		printf("sf_path(%d) = %f\n",i+1,sf_path(i+1));
	}
	sf_param.resize(n_param);
	for (Index i = 1; i <= n_param; i++) {
		if(lb_param(i) != 0 || ub_param(i) != 0) {
			sf_param(i) = max(fabs(lb_param(i)), fabs(ub_param(i)));
		}
		else {
			sf_param(i) = 1;
		}
//		printf("sf_param(%d) = %f\n",i+1,sf_param(i+1));
	}
	sf_events.resize(n_events);
	for (Index i = 1; i <= n_events; i++) {
		if(lb_events(i) != 0 || ub_events(i) != 0) {
			sf_events(i) = max(fabs(lb_events(i)), fabs(ub_events(i)));
		}
		else {
			sf_events(i) = 1;
		}
//		printf("sf_events(%d) = %f\n",i+1,sf_events(i+1));
	}
	sf_t.resize(1);
	if(lb_t0 != 0 || ub_tf != 0) {
		sf_t(1) = max(fabs(lb_t0), fabs(ub_tf));
	}
	else {
		sf_t(1) = 1;
	}
	sf_t(1) = max(fabs(lb_t0), fabs(ub_tf));
}


void OCP::OCPBounds2NLPBounds() {
	determine_scaling_factors();

	Index NLP_n 		= myadolc_nlp->getNLP_n();
	Index NLP_m 		= myadolc_nlp->getNLP_m();
	SVector<double> x_l(NLP_n);
	SVector<double> x_u(NLP_n);
	SVector<double> g_l(NLP_m);
	SVector<double> g_u(NLP_m);
	SVector<double> x_guess(NLP_n);
	SVector<double> node_str(n_nodes - 1);
	SVector<double> sf_NLP_x(NLP_n);
	SVector<double> sf_NLP_g(NLP_m);

	if ((guess.nodes.getRowDim() != (uint)n_nodes) 	||
		(guess.param.getRowDim() != (uint)n_param) 	||
		(guess.u.getColDim() != (uint)n_controls)	||
		(guess.x.getColDim() != (uint)n_states)	||
		(guess.u.getRowDim() != (uint)n_nodes)	||
		(guess.x.getRowDim() != (uint)n_nodes)	) {

		cout<<"===== No guess provided or guess invalid =====\n"
			  "=====  Generating guess based on Bounds  =====\n\n";
		auto_guess_gen();
	}

	// structure of x
	// [..x_1(t_k),..x_n(t_k), u_1(t_k)..u_m(t_k), x_1(t_k+1), ..x_n(t_k+1), u_1(t_k+1), ..u_m(t_k+1), p_1,..p_nparam..]

	Index idx = 1;
	for (Index i = 1; i <= n_nodes; i += 1)
	{
		for (Index j = 1; j <= n_states; j += 1)
		{
			x_l(idx)		= lb_states(j)/sf_x(j);
			x_u(idx)		= ub_states(j)/sf_x(j);
			sf_NLP_x(idx)	= sf_x(j);
			x_guess(idx)	= guess.x(i,j)/sf_x(j);
			idx++;
		}
		for (Index j = 1; j <= n_controls; j += 1)
		{
			x_l(idx)		= lb_controls(j)/sf_u(j);
			x_u(idx)		= ub_controls(j)/sf_u(j);
			sf_NLP_x(idx)	= sf_u(j);
			x_guess(idx)	= guess.u(i,j)/sf_u(j);
			idx++;
		}
	}

	for (Index i = 1; i <= n_param; i += 1)
	{
		x_l(idx)		= lb_param(i)/sf_param(i);
		x_u(idx)		= ub_param(i)/sf_param(i);
		sf_NLP_x(idx)	= sf_param(i);
		x_guess(idx)	= guess.param(i,1)/sf_param(i);
		idx++;
	}

	x_l(idx)		= lb_t0/sf_t(1);
	x_u(idx)		= ub_t0/sf_t(1);
	sf_NLP_x(idx)	= sf_t(1);
	x_guess(idx)	= guess.nodes(1,1)/sf_t(1);
	idx++;

	x_l(idx)		= lb_tf/sf_t(1);
	x_u(idx)		= ub_tf/sf_t(1);
	sf_NLP_x(idx)	= sf_t(1);
	x_guess(idx)	= guess.nodes(n_nodes,1)/sf_t(1);

	double delta_t 	= guess.nodes(n_nodes,1) - guess.nodes(1,1);

	for (Index i = 1; i <= n_nodes - 1; i++) {
		node_str(i)	= (guess.nodes(i+1,1) - guess.nodes(i,1))/delta_t;
	}

  	// structure of g
  	// [..defect(t_k), h1(t_k), .. h_npath(t_k), defect(t_k+1), h1(t_k+1), ..h_npath(t_k+1), events_1, ..events_nevents,..]
	idx = 1;
	for (Index i =1; i <= n_nodes; i += 1)
	{
		for (Index j = 1; j <= n_path; j += 1)
		{
			g_l(idx)		= lb_path(j)/sf_path(j) + 1;
			g_u(idx)		= ub_path(j)/sf_path(j) + 1;
			sf_NLP_g(idx)	= sf_path(j);
			idx++;
		}
		for (Index j = 1; j <= n_states; j += 1)
		{
			if(i < n_nodes) {
				g_l(idx)		= 1.0;
				g_u(idx)		= 1.0;
				sf_NLP_g(idx)	= 1.0;
				idx++;
			}
		}
	}
	for (Index i = 1; i <= n_events; i += 1)
	{
		g_l(idx)	= lb_events(i)/sf_events(i) + 1;
		g_u(idx)	= ub_events(i)/sf_events(i) + 1;
		sf_NLP_g(idx)	= sf_events(i);
		idx++;
	}
/*
	for( Index i = 0; i < NLP_n; i++) {
		sf_NLP_x[i] 	= 1.0;
	}
	for( Index i = 0; i < NLP_m; i++) {
		sf_NLP_g[i]		= 1.0;
	}
	*/
	myadolc_nlp->setBounds(x_l, x_u, g_l, g_u);
	myadolc_nlp->setSF(sf_NLP_x, sf_NLP_g);
	myadolc_nlp->setguess(x_guess);
	myadolc_nlp->setnodestr(node_str);

}

void OCP::auto_guess_gen() {
	SVector<double> nodes = linspace<double>((lb_t0+ub_t0)/2,(lb_tf+ub_tf)/2,n_nodes);
	guess.nodes = nodes;

	SVector<double>* states = new SVector<double>[n_states];
	for (Index i = 0; i < n_states; i++) {
			states[i] = linspace<double>(lb_states(i+1), ub_states(i+1),n_nodes);
		}
	guess.x		= states[0];
	for (Index i = 1; i < n_states; i++) {
		guess.x = (guess.x,states[i]);
	}

	SVector<double>* controls = NULL;
	if (n_controls > 0) {
		controls = new SVector<double>[n_controls];
		for (Index i = 0; i < n_controls; i++) {
			controls[i] = linspace<double>(lb_controls(i+1), ub_controls(i+1),n_nodes);
		}
		guess.u		= controls[0];
		for (Index i = 1; i < n_controls; i++) {
			guess.u = (guess.u,controls[i]);
		}
	}
	else {
		SVector<double> Null_vector;
		guess.u = Null_vector;
	}

	SVector<double> param;
	if(n_param > 0) {
		param.resize((uint) n_param);
		for (Index i = 0; i < n_param; i++) {
			param(i+1) 	= (double)(lb_param(i+1) + ub_param(i+1))/2;
		}
	}
	guess.param = param;

	delete[] states;
	if (controls != NULL) {
		delete[] controls;
	}
}

void OCP::set_endpoint_cost(double (*d_e_cost)(	  const  double* ini_states, const  double* fin_states, const  double* param, const  double& t0, const  double& tf, uint phase),
						   adouble (*ad_e_cost)(const adouble* ini_states, const adouble* fin_states, const adouble* param, const adouble& t0, const adouble& tf, uint phase)) {
	myadolc_nlp->d_e_cost 	= d_e_cost;
	myadolc_nlp->ad_e_cost 	= ad_e_cost;
}
void OCP::set_lagrange_cost( double (*d_l_cost)	( const  double *states, const  double *controls, const  double *param, const  double &time,	uint phase),
							adouble (*ad_l_cost)( const adouble *states, const adouble *controls, const adouble *param, const adouble &time,	uint phase)){
	myadolc_nlp->d_l_cost 	= d_l_cost;
	myadolc_nlp->ad_l_cost 	= ad_l_cost;
}
void OCP::set_derivatives(void (*d_derv)(  double *states_dot,  double *path, const  double *states, const  double *controls, const  double *param, const  double &time, uint phase),
					 void (*ad_derv)(adouble *states_dot, adouble *path, const adouble *states, const adouble *controls, const adouble *param, const adouble &time, uint phase)){
	myadolc_nlp->d_derv 	= d_derv;
	myadolc_nlp->ad_derv 	= ad_derv;

}
void OCP::set_events(void (*d_events)(  double *events, const  double *ini_states, const  double *fin_states, const  double *param, const  double &t0, const  double &tf, uint phase),
				void (*ad_events)(adouble *events, const adouble *ini_states, const adouble *fin_states, const adouble *param, const adouble &t0, const adouble &tf, uint phase)) {
	myadolc_nlp->d_events 	= d_events;
	myadolc_nlp->ad_events  = ad_events;
}
Guess::Guess() {

}

Guess::~Guess() {

}

Config::Config() {
	max_iter 	= 5000;
	NLP_solver 	= ma27;
	warmstart 	= false;
	NLP_tol		= 1e-6;
	opt_oder	= first_order;
	with_mgl	= false;
}
