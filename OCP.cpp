/*
 * OCP.cpp
 *
 *  Created on: May 16, 2015
 *      Author: zineus
 */


#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include "OCP.hpp"

using namespace Ipopt;

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

	lb_states 		= NULL;
	ub_states 		= NULL;
	lb_controls 	= NULL;
	ub_controls 	= NULL;
	lb_param		= NULL;
	ub_param 		= NULL;
	lb_path			= NULL;
	ub_path			= NULL;
	lb_events		= NULL;
	ub_events 		= NULL;

	lb_t0 			= 0;
	ub_t0			= 0;

	lb_tf			= 0;
	ub_tf			= 0;

}

/*
template<class T> void MyADOLC_sparseNLP::euler(const T *states_0, const T *states_dot, const double delta, T *states_1) {
	for (Index i = 0; i < n_states; i += 1)
	{
		states_1[i]		= states_0[i] + delta*states_dot[i];
	}
}

template<class T> void MyADOLC_sparseNLP::trapezoidal(const T *states_0, const T* states_dot_0, const T *states_dot_1, const double delta, T *states_1) {
	for (Index i = 0; i < n_states; i += 1) {
		states_1[i]		= states_0[i] + delta*0.5*(states_0[i] + states_1[i]);
	}
}*/


OCP::~OCP() {
	// TODO Auto-generated destructor stub
	if (lb_states != NULL) {
		delete[] lb_states;
		delete[] ub_states;
		lb_states = NULL;
		ub_states = NULL;
	}
	if (lb_controls != NULL) {
		delete[] lb_controls;
		delete[] ub_controls;
		lb_controls = NULL;
		ub_controls = NULL;
	}
	if (lb_param != NULL) {
		delete[] lb_param;
		delete[] ub_param;
		lb_param = NULL;
		ub_param = NULL;
	}
	if (lb_path != NULL) {
		delete[] lb_path;
		delete[] ub_path;
		lb_path = NULL;
		ub_path = NULL;
	}
	if (lb_events != NULL) {
		delete[] lb_events;
		delete[] ub_events;
		lb_events = NULL;
		ub_events = NULL;
	}
}

ApplicationReturnStatus OCP::set_OCP_structure() {
	Index NLP_n = (n_states + n_controls)*n_nodes + n_param + 2;
//	Index NLP_n = (n_states + n_controls)*n_nodes + n_param;
	Index NLP_m = ((n_nodes - 1)*n_states + n_events + n_path*n_nodes);
	Index* structure = new Index[8];
	structure[0] = n_phases;
	structure[1] = n_nodes;
	structure[2] = n_states;
	structure[3] = n_controls;
	structure[4] = n_param;
	structure[5] = n_events;
	structure[6] = n_path;
	structure[7] = n_linkages;

	myadolc_nlp->setNLP_structure(NLP_n, NLP_m, structure);

	if (n_states > 0) {
		lb_states 	= new Number[n_states];
		ub_states 	= new Number[n_states];
	}

	if (n_controls > 0) {
		lb_controls 	= new Number[n_controls];
		ub_controls 	= new Number[n_controls];
	}

	if (n_param > 0) {
		lb_param 	= new Number[n_param];
		ub_param 	= new Number[n_param];
	}

	if (n_path > 0) {
		lb_path 	= new Number[n_path];
		ub_path 	= new Number[n_path];
	}

	if (n_events > 0) {
		lb_events 		= new Number[n_events];
		ub_events 		= new Number[n_events];
	}

//	app->Options()->SetNumericValue("tol", 1e-6);
//	app->Options()->SetStringValue("mu_strategy", "adaptive");
//	app->Options()->SetStringValue("output_file", "ipopt.out");
//	app->Options()->SetStringValue("nlp_scaling_method","gradient-based");
//	app->Options()->SetStringValue("linear_solver", "ma27");	ma86 & ma57 with memmory leakage
	app->Options()->SetIntegerValue("max_iter", 1000);
	app->Options()->SetStringValue("hessian_approximation", "limited-memory");
  	status = app->Initialize();
  	if (status != Solve_Succeeded) {
  		printf("\n\n*** Error during initialization!\n");
  	}
  	return status;
}

ApplicationReturnStatus OCP::NLP_solve() {

	OCPBounds2NLPBounds();
  	status = app->OptimizeTNLP(myadolc_nlp);

  	if (status == Solve_Succeeded) {
    // Retrieve some statistics about the solve
		Index iter_count = app->Statistics()->IterationCount();
		printf("\n\n*** The problem solved in %d iterations!\n", iter_count);

		Number final_obj = app->Statistics()->FinalObjective();
    	printf("\n\n*** The final value of the objective function is %e.\n", final_obj);
	}
  	Number * NLP_x = myadolc_nlp->getx_opt();
	double t0					= 0.0;
	double tf 					= NLP_x[myadolc_nlp->getNLP_n() - 1];
	double delta 				= (tf - t0)/((double)n_nodes-1);

	nodes_opt.resize(1, n_nodes);
	x_opt.resize(n_nodes, n_states);
	u_opt.resize(n_nodes, n_controls);
	param_opt.resize(1, n_param);

	nodes_opt(1,1) = t0;
	for (Index i = 2; i <= n_nodes; i++) {
		nodes_opt(1,i) 	= nodes_opt(1,i-1) + delta;
	}
	nodes_opt.save("nodes_opt.dat");

	Index idx = 0;
	while (idx < n_nodes*(n_states + n_controls) + n_param) {
		for (Index i = 1; i <= n_nodes; i += 1)	{
			for (Index j = 1; j <= n_states; j += 1){
				x_opt(i,j) 	= NLP_x[idx];
				idx++;
			}
			for (Index j = 1; j <= n_controls; j += 1){
				u_opt(i,j) 	= NLP_x[idx];
				idx++;
			}
		}
		for (Index i = 1; i <= n_param; i += 1)
		{
			param_opt(1,i)			= NLP_x[idx];
			idx++;
		}
	}
	x_opt.save("x_opt.dat");
	u_opt.save("u_opt.dat");
	param_opt.save("param_opt.dat");

	x0_opt.resize(1,n_states);
	xf_opt.resize(1,n_states);
	for (Index i = 1; i <= n_states; i += 1) {
		x0_opt(1,i) 		= x_opt(1,i);
		xf_opt(1,i)			= x_opt(n_nodes,i);
	}
	x0_opt.save("x0_opt.dat");
	xf_opt.save("xf_opt.dat");

  	return status;
}

void OCP::OCPBounds2NLPBounds() {
	Index NLP_n 		= myadolc_nlp->getNLP_n();
	Index NLP_m 		= myadolc_nlp->getNLP_m();
	Number* x_l 		= new Number[NLP_n];
	Number* x_u 		= new Number[NLP_n];
	Number* g_l 		= new Number[NLP_m];
	Number* g_u 		= new Number[NLP_m];
	Number* x_guess 	= new Number[NLP_n];
	Number* node_str 	= new Number[n_nodes - 1];
//	Number* sf_NLP_x	= new Number[NLP_n];

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

	Index idx = 0;
	for (Index i = 0; i < n_nodes; i += 1)
	{
		for (Index j = 0; j < n_states; j += 1)
		{
			x_l[idx]		= lb_states[j];
			x_u[idx]		= ub_states[j];
//			sf_NLP_x[idx]	= get_sf(lb_states[j], ub_states[j]);
			x_guess[idx]	= guess.x(i+1,j+1);
			idx++;
		}
		for (Index j = 0; j < n_controls; j += 1)
		{
			x_l[idx]		= lb_controls[j];
			x_u[idx]		= ub_controls[j];
			x_guess[idx]	= guess.u(i+1,j+1);
			idx++;
		}
	}

	for (Index i = 0; i < n_param; i += 1)
	{
		x_l[idx]		= lb_param[i];
		x_u[idx]		= ub_param[i];
		x_guess[idx]	= guess.param(1,i+1);
		idx++;
	}

	x_l[idx]		= lb_t0;
	x_u[idx]		= ub_t0;
	x_guess[idx]	= guess.nodes(1,1);
	idx++;

	x_l[idx]		= lb_tf;
	x_u[idx]		= ub_tf;
	x_guess[idx]	= guess.nodes(n_nodes,1);

	double delta_t 	= guess.nodes(n_nodes,1) - guess.nodes(1,1);

	for (Index i = 0; i < n_nodes - 1; i++) {
		node_str[i]	= (guess.nodes(i+2,1) - guess.nodes(i+1,1))/delta_t;
	}

  	// structure of g
  	// [..defect(t_k), h1(t_k), .. h_npath(t_k), defect(t_k+1), h1(t_k+1), ..h_npath(t_k+1), events_1, ..events_nevents,..]
	idx = 0;
	for (Index i = 0; i < n_nodes; i += 1)
	{
		for (Index j = 0; j < n_path; j += 1)
		{	cout<<"path idx = "<<idx<<"\n";
			g_l[idx]	= lb_path[j];
			g_u[idx]	= ub_path[j];
			idx++;
		}
		for (Index j = 0; j < n_states; j += 1)
		{
			if(i < n_nodes - 1) {
				g_l[idx]	= 0.0;
				g_u[idx]	= 0.0;
				idx++;
			}
		}
	}
	for (Index i = 0; i < n_events; i += 1)
	{
		g_l[idx]	= lb_events[i];
		g_u[idx]	= ub_events[i];
		idx++;
	}
	myadolc_nlp->setBounds(x_l, x_u, g_l, g_u);
	myadolc_nlp->setguess(x_guess);
	myadolc_nlp->setnodestr(node_str);

}

void OCP::auto_guess_gen() {
	SVector nodes = linspace((lb_t0+ub_t0)/2,(lb_tf+ub_tf)/2,n_nodes);
	guess.nodes = nodes;

	SVector* states = new SVector[n_states];
	for (Index i = 0; i < n_states; i++) {
			states[i] = linspace(lb_states[i], ub_states[i],n_nodes);
		}
	guess.x		= states[0];
	for (Index i = 1; i < n_states; i++) {
		guess.x = (guess.x,states[i]);
	}

	SVector* controls = NULL;
	if (n_controls > 0) {
		controls = new SVector[n_controls];
		for (Index i = 0; i < n_controls; i++) {
			controls[i] = linspace(lb_controls[i], ub_controls[i],n_nodes);
		}
		guess.u		= controls[0];
		for (Index i = 1; i < n_controls; i++) {
			guess.u = (guess.u,controls[i]);
		}
	}
	else {
		SVector Null_vector;
		guess.u = Null_vector;
	}

	SVector param;
	if(n_param > 0) {
		param.resize((uint) n_param);
		for (Index i = 0; i < n_param; i++) {
			param.setelement((uint)i,(double)(lb_param[i] + ub_param[i])/2);
		}
	}
	guess.param = param;

	delete[] states;
	if (controls != NULL) {
		delete[] controls;
	}
}

Guess::Guess() {

}

Guess::~Guess() {
//	delete[] nodes; delete[] x, delete[] u, delete[] param;
}
