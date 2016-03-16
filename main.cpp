#include "ADOL-C_sparseNLP.hpp"


#define N_NODES		15

template<class T> T endpoint_cost (	const T* ini_states,
									const T* fin_states,
									const T* param,
									const T& t0,
									const T& tf,
									Index phase, const double* constants) {


	return fin_states[2];

}


template<class T> void derivatives(	T *states_dot,
									T *path,
							 	 	const T *states,
							 	 	const T *controls,
							 	 	const T *param,
							 	 	const T &time,
							 	 	Index phase, const double* constants) {
//	T x 	= states[0];
	T v 	= states[1];
//	T cost 	= states[2];

	T u = controls[0];

	states_dot[0] 	= v;
	states_dot[1]	= u;
	states_dot[2]	= u*u/2.;
}

template<class T> void events(	T *events,
								const T *ini_states,
								const T *fin_states,
								const T *param,
								const T &t0,
								const T &tf,
								Index phase, const double* constants) {

	events [0]	= ini_states[0];
	events [1]	= ini_states[1];
	events [2]	= fin_states[0];
	events [3]	= fin_states[1];
	events [4]	= ini_states[2];

}

using namespace Ipopt;


int main(int argv, char* argc[])
{
	ApplicationReturnStatus status;
	SmartPtr<MyADOLC_sparseNLP> problem = new MyADOLC_sparseNLP();
	SmartPtr<IpoptApplication>app 		= new IpoptApplication();

	problem->n_nodes 			= N_NODES;
	problem->n_states 			= 3;
	problem->n_controls			= 1;
	problem->n_parameters		= 0;
	problem->n_events 			= 5;
	problem->n_path_constraints	= 0;

	problem->set_endpoint_cost(&endpoint_cost, &endpoint_cost);
	problem->set_events(&events, &events);
	problem->set_derivatives(&derivatives, &derivatives);

	problem->config.max_iter 		= 5000;
	problem->config.NLP_solver 		= ma27;
	problem->config.warmstart 		= false;
	problem->config.NLP_tol			= 1e-6;
	problem->config.with_mgl		= false;
	problem->config.disc_method		= Hermite_Simpson;
	problem->config.H_approximation = false;

	problem->mem_allocation();

	problem->lb_states(1)	= 0.0;
	problem->lb_states(2)	= -10.0;
	problem->lb_states(3)	= 0.0;

	problem->ub_states(1)	= 2.0;
	problem->ub_states(2)	= 10.0;
	problem->ub_states(3)	= 10.0;

	problem->lb_controls(1)	=-4.0;
	problem->ub_controls(1)	= 4.0;

	problem->lb_t0 			= 0.0;
	problem->ub_t0 			= 0.0;

	problem->lb_tf			= 1.0;
	problem->ub_tf			= 1.0;

	problem->lb_events(1)	= 0.0;
	problem->ub_events(1)	= 0.0;

	problem->lb_events(2)	= 1.0;
	problem->ub_events(2)	= 1.0;
	
	problem->lb_events(3)	= 0.0;
	problem->ub_events(3)	= 0.0;

	problem->lb_events(4)	= -1.0;
	problem->ub_events(4)	= -1.0;

	problem->lb_events(5)	= 0.0;
	problem->ub_events(5)	= 0.0;

	status = problem->initialization(app);
	status = problem->solve(app);

	problem->results.x.save("x.dat");
	problem->results.u.save("u.dat");
	problem->results.nodes.save("nodes.dat");

	problem->results.lam_defects.save("lam_defects.dat");
	problem->results.lam_events.save("lam_events.dat");

	problem->results.Hamiltonian.save("Hamiltonian.dat");

	return status;

}
