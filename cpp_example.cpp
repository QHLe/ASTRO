
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
	problem.n_nodes 		= 30;
	problem.n_states 		= 2;
	problem.n_controls		= 1;
	problem.n_param			= 0;
	problem.n_events 		= 4;
	problem.n_path			= 0;

//	problem.lagrange_cost 	= &lagrange_cost;
//	problem.endpoint_cost 	= &endpoint_cost;

	ApplicationReturnStatus status;

	status = problem.NLP_initialization();

	problem.lb_states[0]	= - 2.0;
	problem.lb_states[1]	= -10.0;

	problem.ub_states[0]	=  1.0;
	problem.ub_states[1]	= 10.0;

	problem.lb_controls[0]	= -1;
	problem.ub_controls[0]	=  1;

	problem.lb_tf			= 0.1;
	problem.ub_tf			= 10;

	problem.lb_events[0]	= 0.0;
	problem.ub_events[0]	= 0.0;
	problem.lb_events[1]	= 0.0;
	problem.ub_events[1]	= 0.0;
	problem.lb_events[2]	= 1.0;
	problem.ub_events[2]	= 1.0;
	problem.lb_events[3]	= 0.0;
	problem.ub_events[3]	= 0.0;

	status = problem.NLP_solve();

	return (int) status;

}
