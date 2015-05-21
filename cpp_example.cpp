
#include "OCP.hpp"
#include "ADOL-C_sparseNLP.hpp"
#include "cpp_example.hpp"

using namespace Ipopt;


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
	problem.lb_states[1]	=-10.0;
	problem.lb_states[2]	=-10.0;

	problem.ub_states[0]	= 1.0/9.0;
	problem.ub_states[1]	= 10.0;
	problem.ub_states[2]	= 10.0;

	problem.lb_controls[0]	=-10.0;
	problem.ub_controls[0]	= 10.0;

	problem.lb_t0 			= 0.0;
	problem.ub_t0 			= 0.0;

	problem.lb_tf			= 0.1;
	problem.ub_tf			= 1.0;

	problem.lb_events[0]	= 0.0;
	problem.ub_events[0]	= 0.0;

	problem.lb_events[1]	= 1.0;
	problem.ub_events[1]	= 1.0;

	problem.lb_events[2]	= 0.0;
	problem.ub_events[2]	= 0.0;

	problem.lb_events[3]	= 0.0;
	problem.ub_events[3]	= 0.0;

	problem.lb_events[4]	=-1.0;
	problem.ub_events[4]	=-1.0;

/*
	problem.guess.nodes		= guess_nodes;
	problem.guess.x			= (guess_x,guess_v);
	problem.guess.u			= guess_u;
*/
	status = problem.NLP_solve();

	return (int) status;

}
