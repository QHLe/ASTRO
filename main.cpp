#include "ADOL-C_sparseNLP.hpp"


#define N_NODES		11

template<class T> T endpoint_cost (	const T* ini_states,
									const T* fin_states,
									const T* param,
									const T& t0,
									const T& tf,
									Index phase, const double* constants) {


	return tf;//fin_states[2];

}


template<class T> void derivatives(	T *states_dot,
									T *path,
							 	 	const T *states,
							 	 	const T *controls,
							 	 	const T *param,
							 	 	const T &time,
							 	 	Index phase, const double* constants) {
	T x2 = states[1];
//	T x3 = states[2];

	T u = controls[0];

	states_dot[0] 	= x2;
	states_dot[1]	= u;
	states_dot[2]	= u*u/2;
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
	events [4] 	= ini_states[2];

}

using namespace Ipopt;


int main(int argv, char* argc[])
{
	ApplicationReturnStatus status;
	SmartPtr<MyADOLC_sparseNLP> problem = new MyADOLC_sparseNLP();

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
	problem->config.NLP_tol			= 1e-8;
	problem->config.with_mgl		= false;
	problem->config.disc_method		= trapezoidal;
	problem->config.H_approximation = false;

	problem->mem_allocation();

	problem->lb_states(1)	= 0.0;
	problem->lb_states(2)	=-10.0;
	problem->lb_states(3)	= 0.0;

	problem->ub_states(1)	= 2.0;
	problem->ub_states(2)	= 10.0;
	problem->ub_states(3)	= 10.0;

	problem->lb_controls(1)	=-1.0;
	problem->ub_controls(1)	= 1.0;

	problem->lb_t0 			= 0.0;
	problem->ub_t0 			= 0.0;

	problem->lb_tf			= 1.0;
	problem->ub_tf			= 3.0;

	problem->lb_events(1)	= -0.0;
	problem->ub_events(1)	= -0.0;

	problem->lb_events(2)	= 0.0;
	problem->ub_events(2)	= 0.0;
	
	problem->lb_events(3)	= 1.0;
	problem->ub_events(3)	= 1.0;

	problem->lb_events(4)	= -0.0;
	problem->ub_events(4)	= -0.0;

	problem->lb_events(5)	= -0.0;
	problem->ub_events(5)	= -0.0;

	status = problem->initialization();

//	status = problem->solve();
	/**/
		SmartPtr<IpoptApplication>app 				= new IpoptApplication();

		app->Options()->SetStringValue("mu_strategy", "adaptive");
	//	app->Options()->SetStringValue("output_file", "ipopt.out");
	//	app->Options()->SetStringValue("nlp_scaling_method","gradient-based");
		app->Options()->SetStringValue("linear_solver", "mumps");//	ma86 & ma57 with memmory leakage
		app->Options()->SetIntegerValue("max_iter", problem->config.max_iter);
		app->Options()->SetIntegerValue("print_level", problem->config.print_level);


		if (problem->config.H_approximation){
			app->Options()->SetStringValue("hessian_approximation", "limited-memory");
		}
		else {
			app->Options()->SetStringValue("hessian_approximation", "exact");
		}

		status = app->Initialize();
		app->OptimizeTNLP(problem);
		problem->results.x.save("x.dat");
		problem->results.u.save("u.dat");
		problem->results.nodes.save("nodes.dat");
	//*/
	problem->results.x.save("x.dat");
	problem->results.u.save("u.dat");
	problem->results.nodes.save("nodes.dat");

	return 0;

}
