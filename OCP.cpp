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

	lagrange_cost 	= NULL;
	endpoint_cost 	= NULL;
	dynamics		= NULL;
	events			= NULL;
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

ApplicationReturnStatus OCP::NLP_initialization() {
	Index NLP_n = (n_states + n_controls)*n_nodes + n_param + 1;
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
  	return status;
}

void OCP::OCPBounds2NLPBounds() {
	Index NLP_n = myadolc_nlp->getNLP_n();
	Index NLP_m = myadolc_nlp->getNLP_m();
	Number* x_l = new Number[NLP_n];
	Number* x_u = new Number[NLP_n];
	Number* g_l = new Number[NLP_m];
	Number* g_u = new Number[NLP_m];

	// structure of x
	// [..x_1(t_k),..x_n(t_k), u_1(t_k)..u_m(t_k), x_1(t_k+1), ..x_n(t_k+1), u_1(t_k+1), ..u_m(t_k+1), p_1,..p_nparam..]

	Index idx = 0;
	for (Index i = 0; i < n_nodes; i += 1)
	{
		for (Index j = 0; j < n_states; j += 1)
		{
			x_l[idx]	= lb_states[j];
			x_u[idx]	= ub_states[j];
			idx++;
		}
		for (Index j = 0; j < n_controls; j += 1)
		{
			x_l[idx]	= lb_controls[j];
			x_u[idx]	= ub_controls[j];
			idx++;
		}
	}
	for (Index i = 0; i < n_param; i += 1)
	{
		x_l[idx]	= lb_param[i];
		x_u[idx]	= ub_param[i];
		idx++;
	}

	x_l[idx]		= lb_tf;
	x_u[idx]		= ub_tf;

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

/*
	for (Index i = 0; i < NLP_n; i += 1)
	{
		printf("x_l[%d] = %f\tx_u[%d] = %f\n",i,x_l[i],i,x_u[i]);
	}
	for (Index i = 0; i < NLP_m; i += 1)
	{
		printf("g_l[%d] = %f\tg_u[%d] = %f\n",i,g_l[i],i,g_u[i]);
	}
*/
	myadolc_nlp->setBounds(x_l, x_u, g_l, g_u);
}

void OCP::auto_guess_gen() {
	SVector nodes = linspace(lb_t0,ub_tf,n_nodes);
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
		SVector* controls = new SVector[n_controls];
		for (Index i = 0; i < n_controls; i++) {
			states[i] = linspace(lb_controls[i], ub_controls[i],n_nodes);
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
