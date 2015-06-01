#include "OCP.hpp"
#include "ADOL-C_sparseNLP.hpp"
#include "model.hpp"

template<class T> T lagrange_cost(	const T *states,
									const T *controls,
									const T *param,
									const T &time,
									uint phase) {
	return 0;
}

template<class T> T endpoint_cost (	const T* ini_states,
									const T* fin_states,
									const T* param,
									const T& t0,
									const T& tf,
									uint phase) {

	return tf;//fin_states[2];
}


template<class T> void derivatives(	T *states_dot,
									T *path,
							 	 	const T *states,
							 	 	const T *controls,
							 	 	const T *param,
							 	 	const T &time,
							 	 	uint phase) {

	T x1 = states[0];
	T x2 = states[1];
//	T x3 = states[2];

	T u = controls[0];

	states_dot[0] 	= x2;
	states_dot[1]	= u;
//	states_dot[2]	= u*u/2;
}

template<class T> void events(	T *events,
								const T *ini_states,
								const T *fin_states,
								const T *param,
								const T &t0,
								const T &tf,
								uint phase) {

/*	events [0]	= ini_states[0];
	events [1]	= ini_states[1];
	events [2]	= ini_states[2];
	events [3]	= fin_states[0];
	events [4]	= fin_states[1];*/

	events [0]	= ini_states[0];
	events [1]	= ini_states[1];
	events [2]	= fin_states[0];
	events [3]	= fin_states[1];

}

using namespace Ipopt;


int main(int argv, char* argc[])
{
	OCP problem;
	problem.n_nodes 		= N_NODES;
	problem.n_states 		= 2;
	problem.n_controls		= 1;
	problem.n_param			= 0;
	problem.n_events 		= 4;
	problem.n_path			= 0;

	problem.set_lagrange_cost(&lagrange_cost, &lagrange_cost);
	problem.set_endpoint_cost(&endpoint_cost, &endpoint_cost);
	problem.set_events(&events, &events);
	problem.set_derivatives(&derivatives, &derivatives);

	ApplicationReturnStatus status;

	status = problem.set_OCP_structure();

	problem.lb_states[0]	= 0.0;
	problem.lb_states[1]	=-10.0;
//	problem.lb_states[2]	=-10.0;

//	problem.ub_states[0]	= 1.0/9.0;
	problem.ub_states[0]	= 1.0;
	problem.ub_states[1]	= 10.0;
//	problem.ub_states[2]	= 10.0;

	problem.lb_controls[0]	=-10.0;
	problem.ub_controls[0]	= 10.0;

	problem.lb_t0 			= 0.0;
	problem.ub_t0 			= 0.0;

	problem.lb_tf			= 0.1;
	problem.ub_tf			= 1.0;

	problem.lb_events[0]	= 0.0;
	problem.ub_events[0]	= 0.0;

	problem.lb_events[1]	= 0.0;
	problem.ub_events[1]	= 0.0;

	problem.lb_events[2]	= 1.0;
	problem.ub_events[2]	= 1.0;

	problem.lb_events[3]	= 0.0;
	problem.ub_events[3]	= 0.0;

/*
	problem.lb_events[4]	=-1.0;
	problem.ub_events[4]	=-1.0;
*/
/*
	problem.guess.nodes		= guess_nodes;
	problem.guess.x			= (guess_x,guess_v);
	problem.guess.u			= guess_u;
*/
	status = problem.NLP_solve();

	return (int) status;

}
