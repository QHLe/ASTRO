

#define N_NODES		300


#include "OCP.hpp"
#include "cpp_example.hpp"

using namespace Ipopt;

template<class T> T lagrange_cost ( const double* states,
									const double* controls,
									const double* param,
									const double& time,
									uint phase) {
	return 0;
}
template<class T> T endpoint_cost (	const double* ini_states,
									const double* fin_states,
									const double* param,
									const double& t0,
									const double& tf,
									uint phase) {
	return 0;
}


int main(int argv, char* argc[])
{
	OCP problem;
	problem.n_nodes 		= N_NODES;
	problem.n_states 		= 3;
	problem.n_controls		= 1;
	problem.n_param			= 0;
	problem.n_events 		= 5;
	problem.n_path			= 0;

//	problem.lagrange_cost 	= &lagrange_cost;
//	problem.endpoint_cost 	= &endpoint_cost;

	ApplicationReturnStatus status;

	status = problem.set_OCP_structure();

	problem.lb_states[0]	= 0.0;
	problem.lb_states[1]	= 0.0;
	problem.lb_states[2]	= 0.0;

	problem.ub_states[0]	= 20.0;
	problem.ub_states[1]	= 20.0;
	problem.ub_states[2]	= 20.0;

	problem.lb_controls[0]	= 0.0;
	problem.ub_controls[0]	= 3.14;

	problem.lb_t0 			= 0.0;
	problem.ub_t0 			= 0.0;

	problem.lb_tf			= 0.1;
	problem.ub_tf			= 10;

	problem.lb_events[0]	= 0.0;
	problem.ub_events[0]	= 0.0;
	problem.lb_events[1]	= 0.0;
	problem.ub_events[1]	= 0.0;
	problem.lb_events[2]	= 0.0;
	problem.ub_events[2]	= 0.0;
	problem.lb_events[3]	= 2.0;
	problem.ub_events[3]	= 2.0;
	problem.lb_events[4]	= 2.0;
	problem.ub_events[4]	= 2.0;

	SVector guess_nodes		= linspace(0,2,N_NODES);
	SVector guess_x			= linspace(0,1,N_NODES);
	SVector guess_v			= linspace(0,1,N_NODES);
	SVector guess_u			= linspace(1,-1,N_NODES);
	double xi=0, vi=0, ui=1, delta = 2/((double)(N_NODES - 1));
	guess_u.setelement(1,ui);
	guess_v.setelement(1,vi);
	guess_x.setelement(1,xi);
	for (uint i = 2; i <= N_NODES; i++) {
		if (i < N_NODES/2) {
			ui = 1;
		}
		else {
			ui = -1;
		}
		vi = vi + ui*delta;
		xi = xi + vi*delta + ui/2*delta*delta;
		guess_u.setelement(i,ui);
		guess_v.setelement(i,vi);
		guess_x.setelement(i,xi);
	}
/*
	problem.guess.nodes		= guess_nodes;
	problem.guess.x			= (guess_x,guess_v);
	problem.guess.u			= guess_u;
*/
	status = problem.NLP_solve();

	return (int) status;

}
